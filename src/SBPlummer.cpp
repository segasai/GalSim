/* -*- c++ -*-
 * Copyright (c) 2012-2016 by the GalSim developers team on GitHub
 * https://github.com/GalSim-developers
 *
 * This file is part of GalSim: The modular galaxy image simulation toolkit.
 * https://github.com/GalSim-developers/GalSim
 *
 * GalSim is free software: redistribution and use in source and binary forms,
 * with or without modification, are permitted provided that the following
 * conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions, and the disclaimer given in the accompanying LICENSE
 *    file.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions, and the disclaimer given in the documentation
 *    and/or other materials provided with the distribution.
 */

//#define DEBUGLOGGING

#include "SBPlummer.h"
#include "SBPlummerImpl.h"

// Define this variable to find azimuth (and sometimes radius within a unit disc) of 2d photons by
// drawing a uniform deviate for theta, instead of drawing 2 deviates for a point on the unit
// circle and rejecting corner photons.
// The relative speed of the two methods was tested as part of issue #163, and the results
// are collated in devutils/external/time_photon_shooting.
// The conclusion was that using sin/cos was faster for icpc, but not g++ or clang++.
#ifdef _INTEL_COMPILER
#define USE_COS_SIN
#endif

// Define this use the Newton-Raphson method for solving the radial value in SBPlummer::shoot
// rather than using OneDimensionalDeviate.
// The relative speed of the two methods was tested as part of issue #163, and the results
// are collated in devutils/external/time_photon_shooting.
// The conclusion was that using OneDimensionalDeviate was universally quite a bit faster.
// However, we leave this option here in case someone has an idea for massively speeding up
// the solution that might be faster than the table lookup.
//#define USE_NEWTON_RAPHSON

#ifdef DEBUGLOGGING
#include <fstream>
//std::ostream* dbgout = &std::cout;
//int verbose_level = 2;
#endif

namespace galsim {

    SBPlummer::SBPlummer(double r0, double flux,
                                 const GSParamsPtr& gsparams) :
        SBProfile(new SBPlummerImpl(r0, flux, gsparams)) {}

    SBPlummer::SBPlummer(const SBPlummer& rhs) : SBProfile(rhs) {}

    SBPlummer::~SBPlummer() {}

    double SBPlummer::getScaleRadius() const
    {
        assert(dynamic_cast<const SBPlummerImpl*>(_pimpl.get()));
        return static_cast<const SBPlummerImpl&>(*_pimpl).getScaleRadius();
    }

    std::string SBPlummer::SBPlummerImpl::serialize() const
    {
        std::ostringstream oss(" ");
        oss.precision(std::numeric_limits<double>::digits10 + 4);
        oss << "galsim._galsim.SBPlummer("<<getScaleRadius()<<", "<<getFlux();
        oss << ", galsim.GSParams("<<*gsparams<<"))";
        return oss.str();
    }

    LRUCache<GSParamsPtr, PlummerInfo> SBPlummer::SBPlummerImpl::cache(
        sbp::max_exponential_cache);

    SBPlummer::SBPlummerImpl::SBPlummerImpl(
        double r0, double flux, const GSParamsPtr& gsparams) :
        SBProfileImpl(gsparams),
        _flux(flux), _r0(r0), _r0_sq(_r0*_r0), _inv_r0(1./r0), _inv_r0_sq(_inv_r0*_inv_r0),
        _info(cache.get(this->gsparams.duplicate()))
    {
        // For large k, we clip the result of kValue to 0.
        // We do this when the correct answer is less than kvalue_accuracy.
        // (1+k^2 r0^2)^-1.5 = kvalue_accuracy
        _ksq_max = (std::pow(this->gsparams->kvalue_accuracy,-1./1.5)-1.);

        // For small k, we can use up to quartic in the taylor expansion to avoid the sqrt.
        // This is acceptable when the next term is less than kvalue_accuracy.
        // 35/16 (k^2 r0^2)^3 = kvalue_accuracy
        _ksq_min = std::pow(this->gsparams->kvalue_accuracy * 16./35., 1./3.);

        _flux_over_2pi = _flux / (2. * M_PI);
        _norm = _flux_over_2pi * _inv_r0_sq;

        dbg<<"Plummer:\n";
        dbg<<"_flux = "<<_flux<<std::endl;
        dbg<<"_r0 = "<<_r0<<std::endl;
        dbg<<"_ksq_max = "<<_ksq_max<<std::endl;
        dbg<<"_ksq_min = "<<_ksq_min<<std::endl;
        dbg<<"_norm = "<<_norm<<std::endl;
        dbg<<"maxK() = "<<maxK()<<std::endl;
        dbg<<"stepK() = "<<stepK()<<std::endl;
    }

    double SBPlummer::SBPlummerImpl::maxK() const
    { return _info->maxK() * _inv_r0; }
    double SBPlummer::SBPlummerImpl::stepK() const
    { return _info->stepK() * _inv_r0; }

    double SBPlummer::SBPlummerImpl::xValue(const Position<double>& p) const
    {
        double r = sqrt(p.x * p.x + p.y * p.y);
        return _norm * std::exp(-r * _inv_r0);
    }

    std::complex<double> SBPlummer::SBPlummerImpl::kValue(const Position<double>& k) const
    {
        double ksq = (k.x*k.x + k.y*k.y)*_r0_sq;

        if (ksq < _ksq_min) {
            return _flux*(1. - 1.5*ksq*(1. - 1.25*ksq));
        } else {
            double temp = 1. + ksq;
            return _flux / (temp * sqrt(temp));
            // NB: flux*std::pow(temp,-1.5) is slower.
        }
    }

    void SBPlummer::SBPlummerImpl::fillXValue(tmv::MatrixView<double> val,
                                                      double x0, double dx, int izero,
                                                      double y0, double dy, int jzero) const
    {
        dbg<<"SBPlummer fillXValue\n";
        dbg<<"x = "<<x0<<" + i * "<<dx<<", izero = "<<izero<<std::endl;
        dbg<<"y = "<<y0<<" + j * "<<dy<<", jzero = "<<jzero<<std::endl;
        if (izero != 0 || jzero != 0) {
            xdbg<<"Use Quadrant\n";
            fillXValueQuadrant(val,x0,dx,izero,y0,dy,jzero);
        } else {
            xdbg<<"Non-Quadrant\n";
            assert(val.stepi() == 1);
            const int m = val.colsize();
            const int n = val.rowsize();
            typedef tmv::VIt<double,1,tmv::NonConj> It;

            x0 *= _inv_r0;
            dx *= _inv_r0;
            y0 *= _inv_r0;
            dy *= _inv_r0;

            for (int j=0;j<n;++j,y0+=dy) {
                double x = x0;
                double ysq = y0*y0;
                It valit = val.col(j).begin();
                for (int i=0;i<m;++i,x+=dx)
                    *valit++ = _norm * std::exp(-sqrt(x*x + ysq));
            }
        }
    }

    void SBPlummer::SBPlummerImpl::fillKValue(tmv::MatrixView<std::complex<double> > val,
                                                      double kx0, double dkx, int izero,
                                                      double ky0, double dky, int jzero) const
    {
        dbg<<"SBPlummer fillKValue\n";
        dbg<<"kx = "<<kx0<<" + i * "<<dkx<<", izero = "<<izero<<std::endl;
        dbg<<"ky = "<<ky0<<" + j * "<<dky<<", jzero = "<<jzero<<std::endl;
        if (izero != 0 || jzero != 0) {
            xdbg<<"Use Quadrant\n";
            fillKValueQuadrant(val,kx0,dkx,izero,ky0,dky,jzero);
        } else {
            xdbg<<"Non-Quadrant\n";
            assert(val.stepi() == 1);
            const int m = val.colsize();
            const int n = val.rowsize();
            typedef tmv::VIt<std::complex<double>,1,tmv::NonConj> It;

            kx0 *= _r0;
            dkx *= _r0;
            ky0 *= _r0;
            dky *= _r0;

            for (int j=0;j<n;++j,ky0+=dky) {
                double kx = kx0;
                double kysq = ky0*ky0;
                It valit = val.col(j).begin();
                for (int i=0;i<m;++i,kx+=dkx) {
                    double ksq = kx*kx + kysq;
                    if (ksq > _ksq_max) {
                        *valit++ = 0.;
                    } else if (ksq < _ksq_min) {
                        *valit++ = _flux * (1. - 1.5*ksq*(1. - 1.25*ksq));
                    } else {
                        double temp = 1. + ksq;
                        *valit++ =  _flux/(temp*sqrt(temp));
                    }
                }
            }
        }
    }

    void SBPlummer::SBPlummerImpl::fillXValue(tmv::MatrixView<double> val,
                                                      double x0, double dx, double dxy,
                                                      double y0, double dy, double dyx) const
    {
        dbg<<"SBPlummer fillXValue\n";
        dbg<<"x = "<<x0<<" + i * "<<dx<<" + j * "<<dxy<<std::endl;
        dbg<<"y = "<<y0<<" + i * "<<dyx<<" + j * "<<dy<<std::endl;
        assert(val.stepi() == 1);
        assert(val.canLinearize());
        const int m = val.colsize();
        const int n = val.rowsize();
        typedef tmv::VIt<double,1,tmv::NonConj> It;

        x0 *= _inv_r0;
        dx *= _inv_r0;
        dxy *= _inv_r0;
        y0 *= _inv_r0;
        dy *= _inv_r0;
        dyx *= _inv_r0;

        It valit = val.linearView().begin();
        for (int j=0;j<n;++j,x0+=dxy,y0+=dy) {
            double x = x0;
            double y = y0;
            for (int i=0;i<m;++i,x+=dx,y+=dyx) *valit++ = _norm * std::exp(-sqrt(x*x + y*y));
        }
    }

    void SBPlummer::SBPlummerImpl::fillKValue(tmv::MatrixView<std::complex<double> > val,
                                                      double kx0, double dkx, double dkxy,
                                                      double ky0, double dky, double dkyx) const
    {
        dbg<<"SBPlummer fillKValue\n";
        dbg<<"kx = "<<kx0<<" + i * "<<dkx<<" + j * "<<dkxy<<std::endl;
        dbg<<"ky = "<<ky0<<" + i * "<<dkyx<<" + j * "<<dky<<std::endl;
        assert(val.stepi() == 1);
        assert(val.canLinearize());
        const int m = val.colsize();
        const int n = val.rowsize();
        typedef tmv::VIt<std::complex<double>,1,tmv::NonConj> It;

        kx0 *= _r0;
        dkx *= _r0;
        dkxy *= _r0;
        ky0 *= _r0;
        dky *= _r0;
        dkyx *= _r0;

        It valit = val.linearView().begin();
        for (int j=0;j<n;++j,kx0+=dkxy,ky0+=dky) {
            double kx = kx0;
            double ky = ky0;
            for (int i=0;i<m;++i,kx+=dkx,ky+=dkyx) {
                double ksq = kx*kx + ky*ky;
                if (ksq > _ksq_max) {
                    *valit++ = 0.;
                } else if (ksq < _ksq_min) {
                    *valit++ = _flux * (1. - 1.5*ksq*(1. - 1.25*ksq));
                } else {
                    double temp = 1. + ksq;
                    *valit++ =  _flux/(temp*sqrt(temp));
                }
            }
        }
    }

    // Constructor to initialize Plummer functions for 1D deviate photon shooting
    PlummerInfo::PlummerInfo(const GSParamsPtr& gsparams)
    {
        dbg<<"Start PlummerInfo with gsparams = "<<gsparams.get()<<std::endl;
#ifndef USE_NEWTON_RAPHSON
        // Next, set up the classes for photon shooting
        _radial.reset(new PlummerRadialFunction());
        dbg<<"Made radial"<<std::endl;
        std::vector<double> range(2,0.);
        range[1] = -std::log(gsparams->shoot_accuracy);
        _sampler.reset(new OneDimensionalDeviate( *_radial, range, true, gsparams));
        dbg<<"Made sampler"<<std::endl;
#endif

        // Calculate maxk:
        _maxk = std::pow(gsparams->maxk_threshold, -1./3.);
        dbg<<"maxk = "<<_maxk<<std::endl;

        // Calculate stepk:
        // int( exp(-r) r, r=0..R) = (1 - exp(-R) - Rexp(-R))
        // Fraction excluded is thus (1+R) exp(-R)
        // A fast solution to (1+R)exp(-R) = x:
        // log(1+R) - R = log(x)
        // R = log(1+R) - log(x)
        double logx = std::log(gsparams->folding_threshold);
        double R = -logx;
        for (int i=0; i<3; i++) R = std::log(1.+R) - logx;
        // Make sure it is at least 5 hlr
        // half-light radius = 1.6783469900166605 * r0
        const double hlr = 1.6783469900166605;
        R = std::max(R,gsparams->stepk_minimum_hlr*hlr);
        _stepk = M_PI / R;
        dbg<<"stepk = "<<_stepk<<std::endl;
    }

    // Set maxK to the value where the FT is down to maxk_threshold
    double PlummerInfo::maxK() const
    { return _maxk; }

    // The amount of flux missed in a circle of radius pi/stepk should be at
    // most folding_threshold of the flux.
    double PlummerInfo::stepK() const
    { return _stepk; }

    boost::shared_ptr<PhotonArray> PlummerInfo::shoot(int N, UniformDeviate ud) const
    {
        dbg<<"PlummerInfo shoot: N = "<<N<<std::endl;
        dbg<<"Target flux = 1.0\n";
        assert(_sampler.get());
        boost::shared_ptr<PhotonArray> result = _sampler->shoot(N,ud);
        dbg<<"PlummerInfo Realized flux = "<<result->getTotalFlux()<<std::endl;
        return result;
    }

    boost::shared_ptr<PhotonArray> SBPlummer::SBPlummerImpl::shoot(
        int N, UniformDeviate u) const
    {
        dbg<<"Plummer shoot: N = "<<N<<std::endl;
        dbg<<"Target flux = "<<getFlux()<<std::endl;
#ifdef USE_NEWTON_RAPHSON
        // The cumulative distribution of flux is 1-(1+r)exp(-r).
        // Here is a way to solve for r by an initial guess followed
        // by Newton-Raphson iterations.  Probably not
        // the most efficient thing since there are logs in the iteration.

        // Accuracy to which to solve for (log of) cumulative flux distribution:
        const double Y_TOLERANCE=this->gsparams->shoot_accuracy;

        double fluxPerPhoton = _flux / N;
        boost::shared_ptr<PhotonArray> result(new PhotonArray(N));

        for (int i=0; i<N; i++) {
            double y = u();
            if (y==0.) {
                // In case of infinite radius - just set to origin:
                result->setPhoton(i,0.,0.,fluxPerPhoton);
                continue;
            }
            // Initial guess
            y = -std::log(y);
            double r = y>2. ? y : sqrt(2.*y);
            double dy = y - r + std::log(1.+r);
            while ( std::abs(dy) > Y_TOLERANCE) {
                r = r + (1.+r)*dy/r;
                dy = y - r + std::log(1.+r);
            }
            // Draw another (or multiple) randoms for azimuthal angle
#ifdef USE_COS_SIN
            double theta = 2. * M_PI * u();
            double sint,cost;
            (theta * radians).sincos(sint,cost);
            double rFactor = r * _r0;
            result->setPhoton(i, rFactor * cost, rFactor * sint, fluxPerPhoton);
#else
            double xu, yu, rsq;
            do {
                xu = 2. * u() - 1.;
                yu = 2. * u() - 1.;
                rsq = xu*xu+yu*yu;
            } while (rsq >= 1. || rsq == 0.);
            double rFactor = r * _r0 / std::sqrt(rsq);
            result->setPhoton(i, rFactor * xu, rFactor * yu, fluxPerPhoton);
#endif
        }
#else
        // Get photons from the PlummerInfo structure, rescale flux and size for this instance
        boost::shared_ptr<PhotonArray> result = _info->shoot(N,u);
        result->scaleFlux(_flux_over_2pi);
        result->scaleXY(_r0);
#endif
        dbg<<"Plummer Realized flux = "<<result->getTotalFlux()<<std::endl;
        return result;
    }
}
