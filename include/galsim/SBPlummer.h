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

#ifndef GalSim_SBPlummer_H
#define GalSim_SBPlummer_H
/** 
 * @file SBPlummer.h @brief SBProfile that implements a 2-d exponential profile.
 */

#include "SBProfile.h"

namespace galsim {

    namespace sbp {

        // How many Plummer profiles to save in the cache
        const int max_exponential_cache = 100;

    }

    /** 
     * @brief Plummer Surface Brightness Profile.  
     *
     * Surface brightness profile with I(r) propto exp[-r/r_0] for some scale-length r_0.  This is a
     * special case of the Sersic profile, but is given a separate class since the Fourier transform
     * has closed form and can be generated without lookup tables.
     */
    class SBPlummer : public SBProfile 
    {
    public:
        /** 
         * @brief Constructor - note that `r0` is scale length, NOT half-light radius `re` as in 
         * SBSersic.
         *
         * @param[in] r0       scale length for the profile that scales as `exp[-(r / r0)]`, NOT the
         *                     half-light radius `re`.
         * @param[in] flux     flux.
         * @param[in] gsparams GSParams object storing constants that control the accuracy of image
         *                     operations and rendering, if different from the default.
         */
        SBPlummer(double r0, double flux, const GSParamsPtr& gsparams);

        /// @brief Copy constructor.
        SBPlummer(const SBPlummer& rhs);

        /// @brief Destructor.
        ~SBPlummer();

        /// @brief Returns the scale radius of the Plummer profile.
        double getScaleRadius() const;

    protected:

        class SBPlummerImpl;

    private:
        // op= is undefined
        void operator=(const SBPlummer& rhs);
    };

}

#endif

