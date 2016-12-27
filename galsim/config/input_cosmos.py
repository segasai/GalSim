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
import galsim
import logging

# This file adds input type cosmos_catalog and gsobject typs COSMOSGalaxy.

# The COSMOSCatalog doesn't need anything special other than registration as a valid input type.
# However, we do make a custom Loader so that we can add a logger line with some information about
# the number of objects in the catalog that passed the initial cuts and other basic catalog info.
from .input import RegisterInputType, InputLoader
class COSMOSLoader(InputLoader):
    def setupImage(self, cosmos_cat, config, base, logger):
        if logger: # pragma: no cover
            # Only report as a warning the first time.  After that, use info.
            first = not base.get('_COSMOSLoader_reported_as_warning',False)
            base['_COSMOSLoader_reported_as_warning'] = True
            if first:
                log_level = logging.WARNING
            else:
                log_level = logging.INFO
            if 'input' in base:
                if 'cosmos_catalog' in base['input']:
                    cc = base['input']['cosmos_catalog']
                    if isinstance(cc,list): cc = cc[0]
                    out_str = ''
                    if 'sample' in cc:
                        out_str += '\n  sample = %s'%cc['sample']
                    if 'dir' in cc:
                        out_str += '\n  dir = %s'%cc['dir']
                    if 'file_name' in cc:
                        out_str += '\n  file_name = %s'%cc['file_name']
                    if out_str != '':
                        logger.log(log_level, 'Using user-specified COSMOSCatalog: %s',out_str)
            logger.info("file %d: COSMOS catalog has %d total objects; %d passed initial cuts.",
                        base['file_num'], cosmos_cat.getNTot(), cosmos_cat.getNObjects())
            if 'gal' in base and 'gal_type' in base['gal']:
                if base['gal']['gal_type']=='parametric':
                    logger.log(log_level,"Using parametric galaxies.")
                else:
                    logger.log(log_level,"Using real galaxies.")

RegisterInputType('cosmos_catalog', COSMOSLoader(galsim.COSMOSCatalog))

# The gsobject type coupled to this is COSMOSGalaxy.

def _BuildCOSMOSGalaxy(config, base, ignore, gsparams, logger):
    """@brief Build a COSMOS galaxy using the cosmos_catalog input item.
    """
    cosmos_cat = galsim.config.GetInputObj('cosmos_catalog', config, base, 'COSMOSGalaxy')

    ignore = ignore + ['num']

    # Special: if galaxies are selected based on index, and index is Sequence or Random, and max
    # isn't set, set it to nobjects-1.
    if 'index' in config:
        galsim.config.SetDefaultIndex(config, cosmos_cat.getNObjects())

    kwargs, safe = galsim.config.GetAllParams(config, base,
        req = galsim.COSMOSCatalog.makeGalaxy._req_params,
        opt = galsim.COSMOSCatalog.makeGalaxy._opt_params,
        single = galsim.COSMOSCatalog.makeGalaxy._single_params,
        ignore = ignore)
    if gsparams: kwargs['gsparams'] = galsim.GSParams(**gsparams)

    if 'index' not in kwargs:
        make_kwargs = {}
        make_kwargs['rng'] = galsim.config.check_for_rng(base, logger, 'COSMOSGalaxy')
        kwargs['index'] = cosmos_cat.selectRandomIndices(1, **make_kwargs)

    if 'gal_type' in kwargs and kwargs['gal_type'] == 'real':
        if 'make_kwargs' in locals():
            kwargs['rng'] = make_kwargs.get('rng',None)
        else:
            kwargs['rng'] = galsim.config.check_for_rng(base, logger, 'COSMOSGalaxy')

    # NB. Even though index is officially optional, it will always be present, either because it was
    #     set by a call to selectRandomIndices or because it was specified explicitly.
    index = kwargs['index']
    if index >= cosmos_cat.getNObjects():
        raise IndexError(
            "%s index has gone past the number of entries in the catalog"%index)

    logger.debug('obj %d: COSMOSGalaxy kwargs = %s',base.get('obj_num',0),kwargs)

    kwargs['cosmos_catalog'] = cosmos_cat

    # Use a staticmethod of COSMOSCatalog to avoid pickling the result of makeGalaxy()
    # The RealGalaxy in particular has a large serialization, so it is more efficient to
    # make it in this process, which is what happens here.
    gal = galsim.COSMOSCatalog._makeSingleGalaxy(**kwargs)

    return gal, safe

# Register this as a valid gsobject type
from .gsobject import RegisterObjectType
RegisterObjectType('COSMOSGalaxy', _BuildCOSMOSGalaxy, input_type='cosmos_catalog')
