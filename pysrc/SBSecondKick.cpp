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

#include "galsim/IgnoreWarnings.h"

#define BOOST_NO_CXX11_SMART_PTR
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"

#include "SBSecondKick.h"

namespace bp = boost::python;

namespace galsim {

    struct PySBSecondKick
    {
        static void wrap()
        {
            bp::class_<SBSecondKick,bp::bases<SBProfile> >("SBSecondKick", bp::no_init)
                .def(bp::init<double,double,double,double,double,double,double,
                             boost::shared_ptr<GSParams> >(
                        (bp::arg("lam"), bp::arg("r0"), bp::arg("L0"), bp::arg("D"), bp::arg("obs"),
                         bp::arg("kcrit"), bp::arg("flux"), bp::arg("gsparams")=bp::object())
                ))
                .def(bp::init<const SBSecondKick &>())
                .def("getLam", &SBSecondKick::getLam)
                .def("getR0", &SBSecondKick::getR0)
                .def("getL0", &SBSecondKick::getL0)
                .def("getD", &SBSecondKick::getD)
                .def("getObs", &SBSecondKick::getObs)
                .def("getKCrit", &SBSecondKick::getKCrit)
                .def("structureFunction", &SBSecondKick::structureFunction)
                .def("tau0", &SBSecondKick::tau0)
                .def("PSF", &SBSecondKick::PSF)
                .enable_pickling()
                ;
        }
    };

    void pyExportSBSecondKick()
    {
        PySBSecondKick::wrap();
    }

} // namespace galsim
