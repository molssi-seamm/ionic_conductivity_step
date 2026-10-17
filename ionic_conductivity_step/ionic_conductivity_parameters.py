# -*- coding: utf-8 -*-
"""
Control parameters for the Ionic Conductivity step in a SEAMM flowchart
"""

import logging
import seamm
from seamm_md_analysis.step_support import source_parameters

logger = logging.getLogger(__name__)


class IonicConductivityParameters(seamm.Parameters):
    """
    The control parameters for Ionic Conductivity.

    The MD to analyze is run earlier in the flowchart -- in a Loop for independent
    runs -- or in other jobs (see seamm_md_analysis.step_support.source_parameters).
    """

    parameters = {
        **source_parameters("the ionic conductivity"),
        "frame interval": {
            "default": 1.0,
            "kind": "float",
            "default_units": "ps",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Longest interval between frames:",
            "help_text": "The positions must be sampled at least this often.",
        },
        "fit start": {
            "default": 0.1,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Fit from:",
            "help_text": (
                "Where the fit of the mean square displacements starts, as a "
                "fraction of their length."
            ),
        },
        "fit end": {
            "default": 0.5,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Fit to:",
            "help_text": (
                "Where the fit ends, as a fraction of the length: the late part has "
                "few time origins and is noisy."
            ),
        },
        "Yeh-Hummer correction": {
            "default": "with the viscosity from a Viscosity step",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("with the viscosity from a Viscosity step", "no"),
            "format_string": "s",
            "description": "Yeh-Hummer correction:",
            "help_text": (
                "Correct the self-diffusion coefficients, and so the Nernst-Einstein "
                "conductivity, for the size of the cell, using the viscosity from a "
                "Viscosity step earlier in the flowchart (skipped, with a note, if "
                "there is none)."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("IonicConductivityParameters.__init__")

        super().__init__(
            defaults={**IonicConductivityParameters.parameters, **defaults}, data=data
        )
