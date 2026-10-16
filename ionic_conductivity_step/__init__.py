# -*- coding: utf-8 -*-

"""
ionic_conductivity_step
A SEAMM plug-in for Ionic Conductivity
"""

# Bring up the classes so that they appear to be directly in
# the ionic_conductivity_step package.

from .ionic_conductivity import IonicConductivity  # noqa: F401, E501
from .ionic_conductivity_parameters import IonicConductivityParameters  # noqa: F401
from .ionic_conductivity_step import IonicConductivityStep  # noqa: F401, E501
from .tk_ionic_conductivity import TkIonicConductivity  # noqa: F401, E501

from .metadata import metadata  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
