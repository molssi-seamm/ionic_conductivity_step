"""Non-graphical part of the Ionic Conductivity step in a SEAMM flowchart

The ionic conductivity of an electrolyte from equilibrium MD with fixed charges:
by Einstein-Helfand, from the collective displacement of the ions' charges,
σ = lim ⟨‖P(t) − P(0)‖²⟩ / (6 V kB t) with P = Σ qᵢ rᵢ over the ions; by
Nernst-Einstein, from the ions' self-diffusion; and their ratio, the ionicity.
The analysis is in :mod:`seamm_md_analysis.electric`.

The MD is run before this step -- in a Loop for independent runs -- or read from
files; the step analyzes it.
"""

from collections import Counter
import importlib.resources
import logging
from pathlib import Path
import warnings

import numpy as np
from seamm_md_analysis import (
    Estimate,
    component_check,
    isotropy,
    open_md,
)
from seamm_md_analysis.electric import (
    cell_dipole,
    classify_molecules,
    conductivity_from_msd,
    nernst_einstein,
)
from seamm_md_analysis.transport import (
    fit_slope,
    msd_by_species,
    squared_displacement_sums,
    yeh_hummer_correction,
)
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    print_details,
    check_md_flowchart,
    preceding_viscosity,
    md_inputs,
    read_files,
    report_state,
    size_warnings,
    state_data,
)
from tabulate import tabulate

import ionic_conductivity_step
import molsystem
import seamm
import seamm_util
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
# MDAnalysis warns about optional packages (netCDF4) in the step's output
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)
job = printing.getPrinter()
printer = printing.getPrinter("Ionic Conductivity")

# Add this module's properties to the standard properties
path = importlib.resources.files("ionic_conductivity_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

#: The fewest cation-anion pairs for a useful conductivity
MIN_ION_PAIRS = 10
#: The fewest frames in a block of a single run
MIN_BLOCK_FRAMES = 20
#: Å²/ps in m²/s
A2_PER_PS = 1.0e-8
#: The unit of the self-diffusion coefficients in the output, m²/s
D_UNIT = 1.0e-9
AXES = ("xx", "yy", "zz")


def hill_formula(symbols):
    """The formula of a molecule in Hill order."""
    counts = Counter(symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(e for e in counts if e not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


def formal_charges(configuration):
    """The formal charge of each atom of the structure, or None if it has none.

    For an MD without fixed charges (a machine-learned forcefield): the ions' net
    charges are all the conductivity needs, and where the charge of an ion sits
    within it does not matter for its displacement.
    """
    atoms = configuration.atoms
    if "formal_charge" not in atoms:
        return None
    q = np.array([c or 0 for c in atoms.get_column_data("formal_charge")], dtype=float)
    return q if np.any(q != 0) else None


def count_ions(charges, molecules):
    """The number of cations and anions, from the charges of the atoms.

    Parameters
    ----------
    charges : array_like
        The charge of each atom, e.
    molecules : [[int]]
        The atoms of each molecule.

    Returns
    -------
    (int, int)
        The numbers of cations and anions.
    """
    classes = classify_molecules(charges, molecules)
    z = [int(classes["z"][m]) for m in classes["ions"]]
    return sum(1 for q in z if q > 0), sum(1 for q in z if q < 0)


def pairs_note(n_cations, n_anions):
    """A warning if there are too few ion pairs for useful statistics, or None."""
    pairs = min(n_cations, n_anions)
    if pairs >= MIN_ION_PAIRS:
        return None
    return (
        f"There are only {pairs} cation-anion pairs ({n_cations} cations, "
        f"{n_anions} anions): the conductivity needs at least about {MIN_ION_PAIRS}, "
        "or it and its error come out poorly. Use a larger system."
    )


class IonicConductivity(seamm.Node):
    """
    The non-graphical part of a Ionic Conductivity step in a flowchart.

    Attributes
    ----------
    parameters : IonicConductivityParameters
        The control parameters for Ionic Conductivity.
    """

    def __init__(
        self,
        flowchart=None,
        title="Ionic Conductivity",
        extension=None,
        logger=logger,
    ):
        logger.debug(f"Creating Ionic Conductivity {self}")
        super().__init__(
            flowchart=flowchart,
            title="Ionic Conductivity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = ionic_conductivity_step.metadata
        self.parameters = ionic_conductivity_step.IonicConductivityParameters()
        self._reset()

    @property
    def version(self):
        """The semantic version of this module."""
        return ionic_conductivity_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return ionic_conductivity_step.__git_revision__

    def _reset(self):
        """Clear the data of earlier runs, e.g. from an earlier loop iteration."""
        self.paths = []  # The md.h5md file of each run

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "ionic-conductivity-step"
        parser = seamm_util.getParser()

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )
        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do."""
        if P is None:
            P = self.parameters.values_to_dict()

        what = (
            "the ionic conductivity by Einstein-Helfand and Nernst-Einstein, and "
            "the ionicity"
        )

        source = P["md source"]
        text = f"Calculate {what}, from "
        if source == read_files:
            text += f"the MD in {P['md files']}."
        else:
            text += f"{source} earlier in the flowchart."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Ionic Conductivity step."""
        next_node = super().run(printer)
        self._reset()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))

        _, configuration = self.get_system_configuration()

        self.references.cite(
            raw=self._bibliography["Hansen2013"],
            alias="Hansen2013",
            module="ionic_conductivity_step",
            level=1,
            note="The Einstein-Helfand conductivity and the Nernst-Einstein relation.",
        )

        self.paths = md_inputs(self, P, name="all/position")
        printer.job(
            __(
                f"Analyzing {len(self.paths)} MD run"
                f"{'s' if len(self.paths) > 1 else ''}: "
                + ", ".join(str(p) for p in self.paths),
                indent=self.indent + 4 * " ",
            )
        )
        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        self.analyze(P=P, configuration=configuration)
        return next_node

    # ------------------------------------------------------------------
    # The analysis
    # ------------------------------------------------------------------
    def _read(self, path, configuration, ion_atoms, ion_charges, ions, masses):
        """One run: the collective charge displacement and the centres of mass of
        the ions in each frame, the interval between frames (ps), the mean volume
        (Å³) and temperature (K), and whether the run is isotropic."""
        with open_md(path) as md:
            if not md.periodic:
                raise RuntimeError(f"The MD in {path} is not periodic.")
            times = md.first["particles/all/position/time"]
            unit = times.attrs.get("unit", "ps")
            if isinstance(unit, bytes):
                unit = unit.decode()
            dt = (
                Q_(float(times[1] - times[0]), unit).m_as("ps")
                if times.shape[0] > 1
                else None
            )
            T = float(md.observable("temperature", units="K").values.mean())
            isotropic = isotropy(md)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = md.universe(configuration)
            P = []
            com = []
            V = []
            weights = [masses[m] / masses[m].sum() for m in ions]
            for ts in u.trajectory:
                x = ts.positions
                P.append(cell_dipole(x, ion_atoms, ion_charges))
                com.append([w @ x[m] for m, w in zip(ions, weights)])
                V.append(ts.volume)
        return {
            "P": np.array(P),
            "com": np.array(com),
            "dt": dt,
            "V": float(np.mean(V)),
            "T": T,
            "isotropic": isotropic,
        }

    @staticmethod
    def _sample(Pc, com, dt, V, T, groups, start, end):
        """σ and the self-diffusion coefficients, by component, of one run or
        block, and the MSDs for the graph.

        Returns
        -------
        dict
            ``sigma``: σ_xx, σ_yy, σ_zz (S/m); ``D``: (n_species, 3), m²/s;
            ``msd``: the collective MSD (e² Å²) and ``species msd`` the MSD of each
            species (Å²), summed over x, y and z; ``V``, ``T``.
        """
        n = len(Pc)
        ts = (np.arange(n) * dt).tolist()
        msd = squared_displacement_sums(Pc) / (n - np.arange(n))[:, None]
        # σ_aa = ⟨Δ P_a²⟩ / (2 V kB t), three times the 3D formula for one axis
        sigma = np.array(
            [
                3 * conductivity_from_msd(msd[:, a], dt, V, T, start, end)["sigma"]
                for a in range(3)
            ]
        )
        species_msd, _ = msd_by_species(com, groups)
        D = np.array(
            [
                [
                    fit_slope(m[:, a], ts, start=start, end=end).slope / 2 * A2_PER_PS
                    for a in range(3)
                ]
                for m in species_msd
            ]
        )
        return {
            "sigma": sigma,
            "D": D,
            "msd": msd.sum(axis=1),
            "species msd": [m[:, 3] for m in species_msd],
            "V": V,
            "T": T,
        }

    def _ions(self, configuration, charges, scale=1.0):
        """The ions in the structure, from the charges of the MD, and their formal
        charges. ``scale`` is the factor the charges of the ions were scaled by."""
        molecules = configuration.find_molecules(as_indices=True)
        symbols = configuration.atoms.symbols
        classes = classify_molecules(charges, molecules, scale=scale)
        ions = [np.array(molecules[m]) for m in classes["ions"]]
        names = [hill_formula([symbols[i] for i in m]) for m in ions]
        z = {name: int(classes["z"][m]) for name, m in zip(names, classes["ions"])}
        return ions, names, z

    def analyze(self, indent="", P=None, configuration=None, **kwargs):
        """Analyze, print, plot and store the ionic conductivity."""
        if configuration is None:
            raise RuntimeError("There is no structure to analyze.")

        with open_md(self.paths[0]) as md:
            charges = md.charges("all")
            masses = md.masses("all")
            scale = md.ion_charge_scale
        if charges is None:
            # A machine-learned forcefield has no fixed charges, but the conductivity
            # needs only the ions' net charges: their formal charges
            charges = formal_charges(configuration)
            if charges is None:
                printer.normal(
                    __(
                        "The MD has no fixed charges, as with a machine-learned "
                        "forcefield, and the structure has no formal charges, so "
                        "there are no ions to follow: skipped.",
                        indent=self.indent + 4 * " ",
                    )
                )
                printer.normal("")
                return
            scale = 1.0
            printer.normal(
                __(
                    "The MD has no fixed charges, as with a machine-learned "
                    "forcefield, so the ions are found from, and carry, the formal "
                    "charges of the structure.",
                    indent=self.indent + 4 * " ",
                )
            )
        charges = np.asarray(charges, dtype=float)
        masses = np.asarray(masses, dtype=float)

        ions, names, z = self._ions(configuration, charges, scale)
        if len(ions) == 0:
            printer.normal(
                __(
                    "There are no ions, so there is no ionic conductivity to "
                    "calculate: skipped.",
                    indent=self.indent + 4 * " ",
                )
            )
            printer.normal("")
            return

        species = list(dict.fromkeys(names))
        groups = [[i for i, n in enumerate(names) if n == s] for s in species]
        counts = {s: len(g) for s, g in zip(species, groups)}
        ion_atoms = np.concatenate(ions)
        ion_charges = charges[ion_atoms]
        n_cations = sum(counts[s] for s in species if z[s] > 0)
        n_anions = sum(counts[s] for s in species if z[s] < 0)
        text = "The ions: " + ", ".join(f"{counts[s]} {s} ({z[s]:+d})" for s in species)
        if scale != 1.0:
            text += (
                f". Their charges were scaled by {scale:g} in the MD, so they carry "
                f"{scale:g} times their formal charges, and both conductivities are "
                "those of this scaled-charge model"
            )
        printer.normal(__(text + ".", indent=self.indent + 4 * " "))
        # The charges the ions carry in the MD, for Nernst-Einstein as for the
        # collective displacement of the charges
        q = {s: z[s] * scale for s in species}

        runs = [
            self._read(p, configuration, ion_atoms, ion_charges, ions, masses)
            for p in self.paths
        ]
        nruns = len(runs)
        start, end = P["fit start"], P["fit end"]
        dt = runs[0]["dt"] or 1.0
        n_frames = sum(len(r["P"]) for r in runs)
        isotropic, why = runs[0]["isotropic"]

        # The independent samples: the runs, cut to the same length, or blocks of a
        # single run. With one run, the values are from the whole run and only the
        # errors from the blocks.
        whole = None
        if nruns > 1:
            source = "runs"
            n = min(len(r["P"]) for r in runs)
            samples = [
                self._sample(
                    r["P"][:n], r["com"][:n], dt, r["V"], r["T"], groups, start, end
                )
                for r in runs
            ]
        else:
            source = "blocks"
            run = runs[0]
            n_blocks = max(
                1, min(int(P["number of blocks"]), len(run["P"]) // MIN_BLOCK_FRAMES)
            )
            samples = [
                self._sample(Pc, com, dt, run["V"], run["T"], groups, start, end)
                for Pc, com in zip(
                    np.array_split(run["P"], n_blocks),
                    np.array_split(run["com"], n_blocks),
                )
            ]
            whole = self._sample(
                run["P"], run["com"], dt, run["V"], run["T"], groups, start, end
            )
            if n_blocks == 1:
                samples = [whole]

        # The viscosity for the Yeh-Hummer correction, if any
        notes = []
        record = None
        if P["Yeh-Hummer correction"] != "no":
            record = preceding_viscosity(self)
            if record is None:
                notes.append(
                    "The Yeh-Hummer correction is skipped because no Viscosity step "
                    "ran before this one."
                )

        def derived(s):
            """The isotropic σ and D, σ_NE and the ionicity of a run or block."""
            result = {"sigma": float(s["sigma"].mean())}
            D = {sp: float(d.mean()) for sp, d in zip(species, s["D"])}
            result["D"] = D
            result["sigma_NE"] = nernst_einstein(D, counts, q, s["V"], s["T"])
            used = result["sigma_NE"]
            if record is not None:
                L = s["V"] ** (1 / 3) * 1.0e-10
                result["dD"] = yeh_hummer_correction(s["T"], L, record["value"])
                result["D, corrected"] = {sp: d + result["dD"] for sp, d in D.items()}
                result["sigma_NE, corrected"] = nernst_einstein(
                    result["D, corrected"], counts, q, s["V"], s["T"]
                )
                used = result["sigma_NE, corrected"]
            result["ionicity"] = result["sigma"] / used if used > 0 else float("nan")
            return result

        values = [derived(s) for s in samples]
        best = derived(whole) if whole is not None else None

        def estimate(get, scale=1.0):
            """The estimate over the samples; with one run, the value of the whole
            run and the error from its blocks."""
            if len(values) < 2:
                result = Estimate(get(values[0], samples[0]) * scale, source=source)
            else:
                result = Estimate.from_samples(
                    [get(v, s) * scale for v, s in zip(values, samples)], source=source
                )
            if best is not None:
                result.value = get(best, whole) * scale
            return result

        sigma = estimate(lambda v, s: v["sigma"])
        sigma_ne = estimate(lambda v, s: v["sigma_NE"])
        ionicity = estimate(lambda v, s: v["ionicity"])
        D = {sp: estimate(lambda v, s, sp=sp: v["D"][sp], 1 / D_UNIT) for sp in species}
        if record is not None:
            sigma_ne_c = estimate(lambda v, s: v["sigma_NE, corrected"])
            D_c = {
                sp: estimate(lambda v, s, sp=sp: v["D, corrected"][sp], 1 / D_UNIT)
                for sp in species
            }
            dD = best["dD"] if best is not None else np.mean([v["dD"] for v in values])

        sigma_components = np.array([s["sigma"] for s in samples])
        D_components = {
            sp: np.array([s["D"][i] for s in samples]) / D_UNIT
            for i, sp in enumerate(species)
        }

        # The table
        D_units = "10⁻⁹ m²/s"
        table = {"Property": [], "Value": [], "± (95%)": [], "Units": []}

        def row(label, value, units=""):
            if isinstance(value, Estimate):
                v, e = value.formatted()
            else:
                v, e = ("" if value is None else f"{value:.4g}"), ""
            table["Property"].append(label)
            table["Value"].append(v)
            table["± (95%)"].append(e)
            table["Units"].append(units)

        def component(values, i, whole_value=None):
            if len(values) < 2:
                return float(values[0, i] if whole_value is None else whole_value)
            result = Estimate.from_samples(values[:, i], source=source)
            if whole_value is not None:
                result.value = float(whole_value)
            return result

        row("σ (Einstein-Helfand)", sigma, "S/m")
        if not isotropic:
            for i, axis in enumerate(AXES):
                row(
                    f"σ_{axis}",
                    component(
                        sigma_components,
                        i,
                        None if whole is None else whole["sigma"][i],
                    ),
                    "S/m",
                )
        row("σ_NE (Nernst-Einstein), in the box", sigma_ne, "S/m")
        if record is not None:
            row("σ_NE (Nernst-Einstein), corrected", sigma_ne_c, "S/m")
        row("Ionicity σ/σ_NE", ionicity)
        for i, sp in enumerate(species):
            row(f"D({sp}), in the box", D[sp], D_units)
            if not isotropic:
                for j, axis in enumerate(AXES):
                    row(
                        f"D({sp})_{axis}",
                        component(
                            D_components[sp],
                            j,
                            None if whole is None else whole["D"][i, j] / D_UNIT,
                        ),
                        D_units,
                    )
            if record is not None:
                row(f"D({sp}), corrected", D_c[sp], D_units)
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "decimal", "decimal", "left"),
        )
        length = len(text.splitlines()[0])
        printer.normal(
            __(
                "\n" + "Ionic conductivity".center(length) + "\n" + text + "\n",
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )

        # What the table means
        lines = [
            f"From {n_frames} frames of {nruns} run{'s' if nruns > 1 else ''}, every "
            f"{dt:.3g} ps, for {len(ions)} ions; the mean square displacements are "
            f"fitted from {start:.0%} to {end:.0%} of their length. "
            + sigma.describe_error()
        ]
        if isotropic:
            check = component_check(sigma_components, names=["σ_xx", "σ_yy", "σ_zz"])
            lines.append(
                "σ and D are the averages of their x, y and z components. "
                + check.text.replace("the components", "the components of σ", 1)
            )
            for sp in species:
                check = component_check(
                    D_components[sp], names=[f"D({sp})_{a}" for a in AXES]
                )
                lines.append(
                    check.text.replace(
                        "the components", f"the components of D({sp})", 1
                    )
                )
        else:
            lines.append(
                f"The run is not treated as isotropic ({why}), so the diagonal "
                "components are given; σ and D are the averages of their components."
            )
        lines.append(
            "The collective displacement of the charges is noisy: long runs, or "
            "several, are needed."
        )
        note = pairs_note(n_cations, n_anions)
        if note is not None:
            notes.append(note)
        if record is not None:
            eta = record["value"] * 1000
            eta_text = f"{eta:.4g}"
            if record.get("stderr") is not None:
                eta_text += f" ± {record['stderr'] * 1000:.2g} (standard error)"
            where = f"step {record['id']} ({record['step']})"
            if record.get("method"):
                where += f", {record['method']}"
            notes.append(
                "The self-diffusion coefficients are corrected for the size of the "
                f"cell (Yeh-Hummer): ΔD = {dD / D_UNIT:.3g} × {D_units}, with the "
                f"viscosity η = {eta_text} mPa·s from {where}. σ_NE, corrected, is "
                "from the corrected coefficients, and the ionicity uses it. The "
                "errors of the corrected values do not include that of the viscosity."
            )
            T_eta = record.get("T")
            T_run = float(np.mean([r["T"] for r in runs]))
            if T_eta is not None and abs(T_eta - T_run) > 2.0:
                notes.append(
                    f"The viscosity is for {T_eta:.1f} K but the MD here is at "
                    f"{T_run:.1f} K, so the correction is approximate."
                )
            if not isotropic:
                notes.append(
                    "The Yeh-Hummer correction is for a cubic cell, so it is "
                    "approximate here."
                )
        for line in lines + notes:
            printer.normal(__(line, indent=8 * " "))
        printer.normal("")
        self._print_samples(samples, values, species, source)

        state = report_state(self, printer, self.paths)

        # The graph: the whole run, or the average of the runs
        if whole is not None:
            plot_msd = whole["msd"]
            plot_species = whole["species msd"]
        else:
            plot_msd = np.mean([s["msd"] for s in samples], axis=0)
            plot_species = np.mean([s["species msd"] for s in samples], axis=0)
        plot_ts = np.arange(len(plot_msd)) * dt
        V = float(np.mean([r["V"] for r in runs]))
        T = float(np.mean([r["T"] for r in runs]))
        fit = conductivity_from_msd(plot_msd, dt, V, T, start=start, end=end)
        self._plot(plot_ts, plot_msd, fit, species, plot_species)

        self.model = self._find_model()
        data = {
            "sigma": sigma.value,
            "sigma_NE": sigma_ne.value,
            "D": {sp: D[sp].value * D_UNIT for sp in species},
            "nruns": nruns,
            "n_frames": n_frames,
        }
        if sigma.stderr is not None:
            data["sigma,stderr"] = sigma.stderr
        if sigma_ne.stderr is not None:
            data["sigma_NE,stderr"] = sigma_ne.stderr
        if np.isfinite(ionicity.value):
            data["ionicity"] = ionicity.value
            if ionicity.stderr is not None:
                data["ionicity,stderr"] = ionicity.stderr
        if all(D[sp].stderr is not None for sp in species):
            data["D,stderr"] = {sp: D[sp].stderr * D_UNIT for sp in species}
        if record is not None:
            data["D, corrected"] = {sp: D_c[sp].value * D_UNIT for sp in species}
            data["sigma_NE, corrected"] = sigma_ne_c.value
        if state is not None:
            data.update(state_data(state))
        self.store_results(configuration=configuration, data=data)

    def _print_samples(self, samples, values, species, source):
        """The values of each run or block, in the detailed output."""
        what = "Run" if source == "runs" else "Block"
        table = {what: [], "σ": [], "σ_xx": [], "σ_yy": [], "σ_zz": []}
        for sp in species:
            table[f"D({sp})"] = []
        for i, (s, v) in enumerate(zip(samples, values), start=1):
            table[what].append(i)
            table["σ"].append(f"{v['sigma']:.3f}")
            for axis, value in zip(AXES, s["sigma"]):
                table[f"σ_{axis}"].append(f"{value:.3f}")
            for sp in species:
                table[f"D({sp})"].append(f"{v['D'][sp] / D_UNIT:.3f}")
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("center",) + ("decimal",) * (len(table) - 1),
        )
        length = len(text.splitlines()[0])
        title = f"By {what.lower()}: σ in S/m, D in 10⁻⁹ m²/s".center(length)
        print_details(
            printer,
            str(__("\n" + title + "\n" + text + "\n", indent=8 * " ", wrap=False)),
        )

    def _plot(self, ts, msd, fit, species, species_msds):
        """The collective MSD of the charges with its fit, and the ions' MSDs."""
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Ionic conductivity",
        )
        plot = figure.add_plot("collective")
        x_axis = plot.add_axis("x", label="Time (ps)")
        y_axis = plot.add_axis("y", label="<|P(t)-P(0)|^2> (e^2 Å^2)", anchor=x_axis)
        x_axis.anchor = y_axis
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="collective MSD",
            x=list(ts),
            xlabel="t",
            xunits="ps",
            y=list(msd),
            ylabel="MSD",
            yunits="e^2 Å^2",
            color="#4dbd74",
        )
        plot.add_trace(
            x_axis=x_axis,
            y_axis=y_axis,
            name="fit",
            x=fit["xs"],
            xlabel="t",
            xunits="ps",
            y=fit["ys"],
            ylabel="fit",
            yunits="e^2 Å^2",
            color="black",
        )
        plot2 = figure.add_plot("ions")
        x2 = plot2.add_axis("x", label="Time (ps)")
        y2 = plot2.add_axis("y", label="MSD (Å^2)", anchor=x2)
        x2.anchor = y2
        for name, m in zip(species, species_msds):
            plot2.add_trace(
                x_axis=x2,
                y_axis=y2,
                name=name,
                x=list(ts),
                xlabel="t",
                xunits="ps",
                y=list(m),
                ylabel="MSD",
                yunits="Å^2",
            )
        figure.grid_plots("collective - ions")
        path = Path(self.directory) / "ionic_conductivity.graph"
        figure.dump(path)
        if "html" in self.options and self.options["html"]:
            figure.template = "line.html_template"
            figure.dump(path.with_suffix(".html"))

    def _find_model(self):
        """The model (forcefield or model chemistry) of the MD, for properties."""
        try:
            if self.variable_exists("_model_chemistry"):
                return self.get_variable("_model_chemistry")["level"]
            ff = self.get_variable("_forcefield")
            return ff.current_forcefield
        except Exception:
            return "unknown"

    # ------------------------------------------------------------------
    # What the MD must provide
    # ------------------------------------------------------------------
    def md_requirements(self, P=None):
        """What the analysis needs from the MD: positions every ps or so, in a
        periodic cell, from dynamics that are not perturbed."""
        interval = (
            P["frame interval"] if P is not None else self._static("frame interval")
        )
        if interval is None or isinstance(interval, str):
            interval = Q_(1.0, "ps")
        elif not hasattr(interval, "m_as"):
            interval = Q_(float(interval), "ps")
        return Requirements(
            needs=[
                Need(
                    "all/position",
                    max_time=interval,
                    why="the displacement of the ions",
                )
            ],
            periodic=True,
            dynamics=True,
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed and is
        set up for the analysis."""
        source = self._static("md source")
        if source is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def _structure_charges(self, configuration):
        """The charges of the atoms known before the MD runs, or None.

        These are the forcefield's charges, once the atoms are typed (the column
        ``charges_<forcefield>``), or else the formal charges.
        """
        atoms = configuration.atoms
        keys = [k for k in atoms.keys() if k.startswith("charges_")]
        if len(keys) > 0:
            key = keys[0]
            try:
                ff = self.get_variable("_forcefield")
                if f"charges_{ff.current_forcefield}" in keys:
                    key = f"charges_{ff.current_forcefield}"
            except Exception:
                pass
            return np.array(atoms.get_column_data(key), dtype=float)
        if "formal_charge" in atoms:
            return np.array(atoms.get_column_data("formal_charge"), dtype=float)
        return None

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze.

        Returns
        -------
        [str]
            Warnings about the structure, e.g. that it is small or has too few ions.
        """
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: a conductivity is a property of the bulk."
            )
        result = size_warnings(configuration)
        charges = self._structure_charges(configuration)
        if charges is not None:
            molecules = configuration.find_molecules(as_indices=True)
            n_cations, n_anions = count_ions(charges, molecules)
            if n_cations + n_anions == 0:
                result.append(
                    "There are no ions, so there will be no ionic conductivity to "
                    "calculate."
                )
            else:
                note = pairs_note(n_cations, n_anions)
                if note is not None:
                    result.append(note)
        return [f"'{self.title}': {w}" for w in result]
