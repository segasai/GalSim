# Copyright (c) 2012-2016 by the GalSim developers team on GitHub
# https://github.com/GalSim-developers
#
# This file is part of GalSim: The modular galaxy image simulation toolkit.
# https://github.com/GalSim-developers/GalSim
#
# GalSim is free software: redistribution and use in source and binary forms,
# with or without modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions, and the disclaimer given in the accompanying LICENSE
#    file.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions, and the disclaimer given in the documentation
#    and/or other materials provided with the distribution.
#

import sys
import os
import numpy as np
from galsim_test_helpers import timer
from scipy.special import j0

try:
    import galsim
except ImportError:
    path, filename = os.path.split(__file__)
    sys.path.append(os.path.abspath(os.path.join(path, "..")))
    import galsim

MOCK_INFINITY = 1e11
lam = 500.0
r0 = 0.2
L0 = 30.0
D = 8.36
obs = 0.61
kcrit = 2*np.pi/r0
flux = 1.0

# power spectrum of phase fluctuations
class Fphi(object):
    def __init__(self, r0, L0):
        self.r0 = r0
        self.L0 = L0

    def __call__(self, kappa):
        return 0.033/0.423*self.r0**(-5./3)*(kappa**2 + self.L0**(-2))**(-11./6)

class Dphi(object):
    def __init__(self, Fphi_, kcrit):
        self.Fphi_ = Fphi_
        self.kcrit = kcrit

    def _integrand(self, rho):
        return lambda kappa: self.Fphi_(kappa)*(1-j0(rho*kappa))*kappa

    def __call__(self, rho):
        if isinstance(rho, np.ndarray):
            return np.array([self(r) for r in rho])

        return 8*np.pi**2 * galsim.integ.int1d(self._integrand(rho), self.kcrit, MOCK_INFINITY)

class Tau0(object):
    def __init__(self, D):
        self.D = D

    def __call__(self, rho):
        roD = rho/self.D
        return 2./np.pi*(np.arccos(roD)-roD*np.sqrt(1-roD**2))

prof = galsim._galsim.SBSecondKick(lam, r0, L0, D, obs, kcrit, flux, None)
Dphi_ = Dphi(Fphi(r0, L0), kcrit)
Tau0_ = Tau0(D)


@timer
def test_structure_function():
    rhos = np.logspace(-2, 1, 20)
    np.testing.assert_almost_equal(
            [prof.structureFunction(r) for r in rhos],
            Dphi_(rhos),
            decimal=3)

@timer
def test_tau0():
    rhos = np.logspace(-2, np.log10(D), 20)
    np.testing.assert_almost_equal(
            [prof.tau0(r) for r in rhos],
            Tau0_(rhos),
            decimal=3)

@timer
def test_PSF():
    prof1 = galsim._galsim.SBSecondKick(lam, r0, 10.0, D, obs, 2*np.pi/r0, flux, None)
    prof2 = galsim._galsim.SBSecondKick(lam, r0, 10.0, D, obs, 2*np.pi/(1.5*r0), flux, None)
    prof3 = galsim._galsim.SBSecondKick(lam, r0, 10.0, D, obs, 0.0, flux, None)
    prof4 = galsim.Airy(lam=lam, diam=D)

    alpha = [0.0, 0.01, 0.02, 0.03]
    psf1 = np.array([prof1.PSF(a/206265) for a in alpha])
    psf2 = np.array([prof2.PSF(a/206265) for a in alpha])
    psf3 = np.array([prof3.PSF(a/206265) for a in alpha])
    psf4 = np.array([prof4.xValue(0,a) for a in alpha])

    print(psf1/psf1[0])
    print(psf2/psf2[0])
    print(psf3/psf3[0])
    print(psf4/psf4[0])


if __name__ == "__main__":
    test_structure_function()
    test_tau0()
    test_PSF()
