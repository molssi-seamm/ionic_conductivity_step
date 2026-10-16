================================
SEAMM Ionic Conductivity Plug-in
================================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/ionic_conductivity_step
   :target: https://github.com/molssi-seamm/ionic_conductivity_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/ionic_conductivity_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/ionic_conductivity_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/ionic_conductivity_step/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/ionic_conductivity_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/ionic_conductivity_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/ionic_conductivity_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/ionic_conductivity_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/ionic_conductivity_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/ionic_conductivity_step.svg
   :target: https://pypi.python.org/pypi/ionic_conductivity_step
   :alt: PyPi VERSION

A SEAMM plug-in for the ionic conductivity of electrolytes from MD.

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/ionic_conductivity_step/index.html
* Code: https://github.com/molssi-seamm/ionic_conductivity_step

Features
--------

* The ionic conductivity of electrolytes by Einstein-Helfand (with the correlations
  between ions) and by Nernst-Einstein, the ionicity, and each ion's self-diffusion,
  with the Yeh-Hummer correction from a Viscosity step.
* Runs its own MD or analyzes the MD earlier in the flowchart, from SEAMM's H5MD files.

Acknowledgements
----------------

This package was created with the `molssi-seamm/cookiecutter-seamm-plugin`_ tool, which
is based on the excellent Cookiecutter_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award CHE-2136142.

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
