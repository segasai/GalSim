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

#ifndef GalSim_SBSecondKick_H
#define GalSim_SBSecondKick_H
/**
 * @file SBSecondKick.h @brief SBProfile for the geometric optics approximation second kick.
 */

#include "SBProfile.h"

namespace galsim {
    class SBSecondKick : public SBProfile
    {
    public:
        /**
         * @brief Constructor.
         *
         * @param[in] lam          Wavelength in nm.
         * @param[in] r0           Fried parameter in m (at given wavelength lam).
         * @param[in] L0           Outer scale in m.
         * @param[in] D            Telescope diameter in m.
         * @param[in] obs          Linear fraction obscuration of telescope pupil.
         * @param[in] kcrit        Critical Fourier mode scale.
         * @param[in] flux         Flux.
         * @param[in] gsparams     GSParams.
         */
         SBSecondKick(double lam, double r0, double L0, double D, double obs, double kcrit,
                      double flux, const GSParamsPtr& gsparams);

        /// @brief Copy constructor
        SBSecondKick(const SBSecondKick& rhs);

        /// @brief Destructor.
        ~SBSecondKick();

        /// Getters
        double getLam() const;
        double getR0() const;
        double getL0() const;
        double getD() const;
        double getObs() const;
        double getKCrit() const;

        double structureFunction(double) const;
        double tau0(double) const;
        double PSF(double) const;

        friend class PSFIntegrand;

    protected:

        class SBSecondKickImpl;

    private:
        // op= is undefined
        void operator=(const SBSecondKick& rhs);
    };
}

#endif
