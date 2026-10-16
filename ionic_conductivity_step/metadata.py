# -*- coding: utf-8 -*-

"""This file contains metadata describing the results from IonicConductivity"""

metadata = {}

"""Description of the computational models for IonicConductivity.

Hamiltonians, approximations, and basis set or parameterizations,
only if appropriate for this code. For example::

    metadata["computational models"] = {
        "Hartree-Fock": {
            "models": {
                "PM7": {
                    "parameterizations": {
                        "PM7": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": True,
                            "code": "mopac",
                        },
                        "PM7-TS": {
                            "elements": "1-60,62-83",
                            "periodic": True,
                            "reactions": True,
                            "optimization": False,
                            "code": "mopac",
                        },
                    },
                },
            },
        },
    }
"""
# metadata["computational models"] = {
# }

"""Description of the IonicConductivity keywords.

(Only needed if this code uses keywords)

Fields
------
description : str
    A human readable description of the keyword.
takes values : int (optional)
    Number of values the keyword takes. If missing the keyword takes no values.
default : str (optional)
    The default value(s) if the keyword takes values.
format : str (optional)
    How the keyword is formatted in the MOPAC input.

For example::
    metadata["keywords"] = {
        "0SCF": {
            "description": "Read in data, then stop",
        },
        "ALT_A": {
            "description": "In PDB files with alternative atoms, select atoms A",
            "takes values": 1,
            "default": "A",
            "format": "{}={}",
        },
    }
"""
# metadata["keywords"] = {
# }

"""Properties that IonicConductivity produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#Psi4#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    "sigma": {
        "description": "The ionic conductivity (Einstein-Helfand)",
        "dimensionality": "scalar",
        "property": "ionic conductivity#Einstein-Helfand#{model}",
        "type": "float",
        "units": "S/m",
    },
    "sigma,stderr": {
        "description": "The standard error of the ionic conductivity",
        "dimensionality": "scalar",
        "property": "ionic conductivity,stderr#Einstein-Helfand#{model}",
        "type": "float",
        "units": "S/m",
    },
    "sigma_NE": {
        "description": "The Nernst-Einstein ionic conductivity",
        "dimensionality": "scalar",
        "property": "ionic conductivity#Nernst-Einstein#{model}",
        "type": "float",
        "units": "S/m",
    },
    "ionicity": {
        "description": "The ionicity, the ratio of the conductivities",
        "dimensionality": "scalar",
        "property": "ionicity#Ionic Conductivity#{model}",
        "type": "float",
        "units": "",
    },
    "D": {
        "description": "The self-diffusion coefficient of each ion",
        "dimensionality": "{species: value}",
        "type": "json",
        "units": "m^2/s",
    },
    "D, corrected": {
        "description": (
            "The self-diffusion coefficient of each ion with the Yeh-Hummer "
            "correction"
        ),
        "dimensionality": "{species: value}",
        "type": "json",
        "units": "m^2/s",
    },
    "T": {
        "description": "The mean temperature",
        "dimensionality": "scalar",
        "type": "float",
        "units": "K",
    },
    "V": {
        "description": "The mean volume",
        "dimensionality": "scalar",
        "type": "float",
        "units": "Å^3",
    },
    "nruns": {
        "description": "The number of MD runs analyzed",
        "dimensionality": "scalar",
        "type": "integer",
    },
    "n_frames": {
        "description": "The number of frames analyzed",
        "dimensionality": "scalar",
        "type": "integer",
    },
}
