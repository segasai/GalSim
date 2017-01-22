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

#ifndef GalSim_SBSecondKickImpl_H
#define GalSim_SBSecondKickImpl_H

#include "SBProfileImpl.h"
#include "SBSecondKick.h"
#include "LRUCache.h"
#include "OneDimensionalDeviate.h"
#include "Table.h"

namespace galsim {
    class SBSecondKick::SBSecondKickImpl : public SBProfileImpl
    {
    public:
        SBSecondKickImpl(double lam, double r0, double L0, double D, double obs, double kcrit,
                         double flux, const GSParamsPtr& gsparams);
        ~SBSecondKickImpl() {}

        // double xValue(const Position<double>& p) const;
        // std::complex<double> kValue(const Position<double>& k) const;

        bool isAxisymmetric() const { return true; }
        bool hasHardEdges() const { return false; }
        bool isAnalyticX() const { return false; }
        bool isAnalyticK() const { return false; }

        double maxK() const;
        double stepK() const;

        Position<double> centroid() const
        { return Position<double>(0., 0.); }

        double getFlux() const { return _flux; }
        double getLam() const {return _lam; }
        double getR0() const {return _r0; }
        double getL0() const {return _L0; }
        double getD() const {return _D; }
        double getObs() const {return _obs; }
        double getKCrit() const {return _kcrit; }
        // double maxSB();// const { return _xnorm * _info->xValue(0.); }
        double maxSB() const { return 1.0; }  // no idea how right/wrong this is.

        /**
         * @brief SBSecondKick photon-shooting is done numerically with `OneDimensionalDeviate`
         * class.
         *
         * @param[in] N Total number of photons to produce.
         * @param[in] ud UniformDeviate that will be used to draw photons from distribution.
         * @returns PhotonArray containing all the photons' info.
         */
        boost::shared_ptr<PhotonArray> shoot(int N, UniformDeviate ud) const;

        double xValue(const Position<double>& p) const
        { throw SBError("SBSecondKick::xValue() is not implemented"); }
        std::complex<double> kValue(const Position<double>& p) const
        { throw SBError("SBSecondKick::kValue() is not implemented"); }

        double structureFunction(double rho) const;
        double tau0(double rho) const;
        double PSF(double alpha) const;

        std::string serialize() const;

    private:

        double _lam;
        double _r0;
        double _L0;
        double _D;
        double _obs;
        double _kcrit;
        double _flux;

        // Copy constructor and op= are undefined.
        SBSecondKickImpl(const SBSecondKickImpl& rhs);
        void operator=(const SBSecondKickImpl& rhs);

        mutable boost::shared_ptr<FluxDensity> _radial;
        mutable boost::shared_ptr<OneDimensionalDeviate> _sampler;
        mutable Table<double,double> _structure_fn;
        mutable Table<double,double> _PSF;

        void _buildStructureFunctionLUT();
    };
}

#endif
