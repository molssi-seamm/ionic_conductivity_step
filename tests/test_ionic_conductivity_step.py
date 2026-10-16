#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ionic_conductivity_step` package."""

import pytest  # noqa: F401
import ionic_conductivity_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = ionic_conductivity_step.IonicConductivity()
    assert str(type(result)) == "<class 'ionic_conductivity_step.ionic_conductivity.IonicConductivity'>"
