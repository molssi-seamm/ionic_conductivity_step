#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `ionic_conductivity_step` package."""

import contextlib

import numpy as np
import pytest  # noqa: F401

import ionic_conductivity_step  # noqa: F401
from ionic_conductivity_step import ionic_conductivity as ic

#: SPC/E-like charges of a water: O, H, H
WATER_Q = [-0.8476, 0.4238, 0.4238]


@pytest.fixture()
def configuration():
    """A periodic box of waters, each 3 Å apart, as a function of the numbers of
    waters, Na+ and Cl-, with the forcefield charges in ``charges_test``."""
    from molsystem.system_db import SystemDB

    db = SystemDB(filename="file:ic_db?mode=memory&cache=shared")

    def make(n_water, n_na=0, n_cl=0):
        system = db.create_system(name="test")
        configuration = system.create_configuration(name="test")
        configuration.periodicity = 3
        configuration.coordinate_system = "Cartesian"
        configuration.cell.parameters = (30.0, 30.0, 30.0, 90, 90, 90)
        x, y, z, atno, q = [], [], [], [], []
        bonds_i, bonds_j = [], []
        n = 0
        for k in range(n_water + n_na + n_cl):
            x0 = 3.0 * (k % 9) + 1
            y0 = 3.0 * ((k // 9) % 9) + 1
            z0 = 3.0 * (k // 81) + 1
            if k < n_water:
                x += [x0, x0 + 0.8, x0 - 0.8]
                y += [y0, y0 + 0.6, y0 + 0.6]
                z += [z0, z0, z0]
                atno += [8, 1, 1]
                q += WATER_Q
                bonds_i += [n, n]
                bonds_j += [n + 1, n + 2]
                n += 3
            else:
                x.append(x0)
                y.append(y0)
                z.append(z0)
                is_na = k < n_water + n_na
                atno.append(11 if is_na else 17)
                q.append(1.0 if is_na else -1.0)
                n += 1
        ids = configuration.atoms.append(x=x, y=y, z=z, atno=atno)
        if len(bonds_i) > 0:
            configuration.bonds.append(
                i=[ids[i] for i in bonds_i], j=[ids[j] for j in bonds_j]
            )
        if "charges_test" not in configuration.atoms:
            configuration.atoms.add_attribute("charges_test", coltype="float")
        configuration.atoms.get_column("charges_test")[:] = q
        return configuration, np.array(q)

    yield make
    db.close()


def test_construction():
    """Just create an object and test its type."""
    result = ionic_conductivity_step.IonicConductivity()
    assert (
        str(type(result))
        == "<class 'ionic_conductivity_step.ionic_conductivity.IonicConductivity'>"
    )


def test_no_own_md():
    """The step has no sub-flowchart and no parameters for running its own MD."""
    node = ionic_conductivity_step.IonicConductivity()
    assert not hasattr(node, "subflowchart")
    for key in ("nruns", "errors"):
        assert key not in node.parameters
    assert node.md_requirements().dynamics


def test_no_ions_skipped(configuration, monkeypatch):
    """Pure water is reported as having no ions and skipped, not an error."""
    config, charges = configuration(4)

    class FakeMD:
        def charges(self, name):
            return charges

        def masses(self, name):
            return np.ones(charges.size)

    @contextlib.contextmanager
    def fake_open_md(path):
        yield FakeMD()

    monkeypatch.setattr(ic, "open_md", fake_open_md)
    stored = []
    node = ionic_conductivity_step.IonicConductivity()
    monkeypatch.setattr(node, "store_results", lambda **kw: stored.append(kw))
    monkeypatch.setattr(node, "_read", lambda *a, **kw: pytest.fail("read the MD"))
    printed = []
    monkeypatch.setattr(ic.printer, "normal", lambda text: printed.append(str(text)))
    node._id = ("1",)
    node.paths = ["md.h5md"]
    P = node.parameters.values_to_dict()
    node.analyze(P=P, configuration=config)
    assert stored == []
    assert any("There are no ions" in t and "skipped" in t for t in printed)


def test_check_configuration_ion_pairs(configuration):
    """Too few ion pairs give a warning before the MD; enough do not."""
    node = ionic_conductivity_step.IonicConductivity()

    config, _ = configuration(20, 3, 3)
    warnings = node.check_configuration(config)
    assert any("only 3 cation-anion pairs" in w for w in warnings)

    config, _ = configuration(20, 12, 12)
    warnings = node.check_configuration(config)
    assert not any("cation-anion pairs" in w for w in warnings)

    config, _ = configuration(20)
    warnings = node.check_configuration(config)
    assert any("no ions" in w for w in warnings)


def test_check_configuration_not_periodic(configuration):
    """A non-periodic system is an error before the MD."""
    node = ionic_conductivity_step.IonicConductivity()
    config, _ = configuration(4, 1, 1)
    config.periodicity = 0
    with pytest.raises(RuntimeError, match="not periodic"):
        node.check_configuration(config)


def test_count_ions():
    molecules = [[0, 1, 2], [3], [4], [5]]
    charges = WATER_Q + [1.0, 1.0, -2.0]
    assert ic.count_ions(charges, molecules) == (2, 1)
    assert ic.pairs_note(2, 1).startswith("There are only 1 cation-anion pairs")
    assert ic.pairs_note(10, 12) is None
