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

// #define DEBUGLOGGING

#include "galsim/IgnoreWarnings.h"

#define BOOST_NO_CXX11_SMART_PTR
#include <boost/math/special_functions/bessel.hpp>

#include "SBSecondKick.h"
#include "SBSecondKickImpl.h"

namespace galsim {
    SBSecondKick::SBSecondKick(double lam, double r0, double L0, double D, double obs, double kcrit,
                               double flux, const GSParamsPtr& gsparams) :
        SBProfile(new SBSecondKickImpl(lam, r0, L0, D, obs, kcrit, flux, gsparams)) {}

    SBSecondKick::SBSecondKick(const SBSecondKick &rhs) : SBProfile(rhs) {}

    SBSecondKick::~SBSecondKick() {}

    double SBSecondKick::getLam() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getLam();
    }

    double SBSecondKick::getR0() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getR0();
    }

    double SBSecondKick::getL0() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getL0();
    }

    double SBSecondKick::getD() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getD();
    }

    double SBSecondKick::getObs() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getObs();
    }

    double SBSecondKick::getKCrit() const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).getKCrit();
    }

    double SBSecondKick::structureFunction(double rho) const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).structureFunction(rho);
    }

    double SBSecondKick::tau0(double rho) const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).tau0(rho);
    }

    double SBSecondKick::PSF(double alpha) const
    {
        assert(dynamic_cast<const SBSecondKickImpl*>(_pimpl.get()));
        return static_cast<const SBSecondKickImpl&>(*_pimpl).PSF(alpha);
    }

    SBSecondKick::SBSecondKickImpl::SBSecondKickImpl(double lam, double r0, double L0, double D,
                                                     double obs, double kcrit, double flux,
                                                     const GSParamsPtr& gsparams) :
        SBProfileImpl(gsparams),
        _lam(lam*1e-9), //nm->m
        _r0(r0),
        _L0(L0),
        _D(D),
        _obs(obs),
        _kcrit(kcrit),
        _flux(flux),
        _structure_fn(Table<double,double>::spline),
        _PSF(Table<double,double>::spline)
    {
        dbg<<"SBSecondKickImpl constructor: gsparams = "<<gsparams.get()<<std::endl;
        dbg<<"this->gsparams = "<<this->gsparams.get()<<std::endl;

        // Make structure function lookup table.
        _buildStructureFunctionLUT();
    }

    double SBSecondKick::SBSecondKickImpl::maxK() const
    { return 2.*M_PI*_D/_lam/206265; /*inverse arcsec*/ }

    double SBSecondKick::SBSecondKickImpl::stepK() const
    { return M_PI/3.0; /*Hack for now.*/}

    std::string SBSecondKick::SBSecondKickImpl::serialize() const
    {
        std::ostringstream oss(" ");
        oss.precision(std::numeric_limits<double>::digits10 + 4);
        oss << "galsim._galsim.SBSecondKick("
            <<getLam()<<", "
            <<getR0()<<", "
            <<getL0()<<", "
            <<getD()<<", "
            <<getObs()<<", "
            <<getKCrit()<<", "
            <<getFlux()<<", galsim.GSParams("<<*gsparams<<"))";
        return oss.str();
    }

    class StructureFunctionIntegrand : public std::unary_function<double,double>
    {
    public:
        StructureFunctionIntegrand(double r0, double L0, double rho) :
            _r0fac(std::pow(r0,-5./3)), _L0sqrinv(1./L0/L0), _rho(rho) {}
        double operator()(double kappa) const {
            return _phasePower(kappa)*(1.0-j0(_rho*kappa))*kappa;
        }
    private:
        double _phasePower(double kappa) const {
            return 0.033/0.423*_r0fac*std::pow(kappa*kappa + _L0sqrinv, -11./6);
        }
        double _r0fac;  // r0^(-5./3)
        double _L0sqrinv;
        double _rho;
    };

    void SBSecondKick::SBSecondKickImpl::_buildStructureFunctionLUT() {
        double dlogrho = gsparams->table_spacing * sqrt(sqrt(gsparams->kvalue_accuracy / 10.));
        dbg<<"Using dlogrho = "<<dlogrho<<std::endl;
        double rhomin = 1e-6;
        double rhomax = 10.0;
        for (double logrho = std::log(rhomin)-0.001; logrho < std::log(rhomax); logrho += dlogrho){
            double rho = std::exp(logrho);
            StructureFunctionIntegrand I(_r0, _L0, rho);
            integ::IntRegion<double> reg(_kcrit, integ::MOCK_INF);
            double val = integ::int1d(I, reg,
                                      gsparams->integration_relerr,
                                      gsparams->integration_abserr);
            val *= 8.0*M_PI*M_PI;
            dbg<<"logrho = "<<logrho<<", I("<<rho<<") = "<<val<<std::endl;
            _structure_fn.addEntry(logrho,val);
        }
    }

    double SBSecondKick::SBSecondKickImpl::structureFunction(double rho) const {
        double logrho = std::log(rho);
        if (logrho > _structure_fn.argMax()) {
            return _structure_fn(_structure_fn.argMax());
        } else return _structure_fn(logrho);
    }

    double SBSecondKick::SBSecondKickImpl::tau0(double rho) const {
        double roD = rho/_D;
        return 2/M_PI*(std::acos(roD) - roD*std::sqrt(1.0 - roD*roD));
    }

    class PSFIntegrand : public std::unary_function<double,double>
    {
    public:
        PSFIntegrand(double alpha, double lam, const SBSecondKick::SBSecondKickImpl* const pimpl) :
            _alpha_o_lam(alpha/lam), _pimpl(pimpl) {}
        double operator()(double rho) const {
            return rho * j0(2*M_PI*_alpha_o_lam*rho) * _pimpl->tau0(rho) * std::exp(-0.5*_pimpl->structureFunction(rho));
        }
    private:
        double _alpha_o_lam;
        const SBSecondKick::SBSecondKickImpl* const _pimpl;
    };

    double SBSecondKick::SBSecondKickImpl::PSF(double alpha) const {
        PSFIntegrand I(alpha, _lam, this);
        integ::IntRegion<double> reg(0.0, _D);
        double val = integ::int1d(I, reg,
                                  gsparams->integration_relerr,
                                  gsparams->integration_abserr);
        return val;
    }

    class SecondKickRadialFunction : public FluxDensity
    {
    public:
        SecondKickRadialFunction() {}
        double operator()(double r) const { return 0.; }
    private:
    };

    boost::shared_ptr<PhotonArray> SBSecondKick::SBSecondKickImpl::shoot(int N, UniformDeviate ud) const
    {
        dbg<<"SecondKick::shoot: N = "<<N<<std::endl;
        dbg<<"Target flux = 1.0"<<std::endl;

        if (!_sampler) {
            std::vector<double> range(2, 0.);
            range[1] = 1.0; // No idea how accurate this is.
            _radial.reset(new SecondKickRadialFunction());
            _sampler.reset(new OneDimensionalDeviate(*_radial, range, true, gsparams));
        }

        assert(_sampler.get());
        boost::shared_ptr<PhotonArray> result = _sampler->shoot(N,ud);
        dbg<<"SecondKick Realized flux = "<<result->getTotalFlux()<<std::endl;
        return result;
    }
}
