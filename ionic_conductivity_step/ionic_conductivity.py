"""Non-graphical part of the Ionic Conductivity step in a SEAMM flowchart

The ionic conductivity of an electrolyte from equilibrium MD with fixed charges:
by Einstein-Helfand, from the collective displacement of the ions' charges,
σ = lim ⟨‖P(t) − P(0)‖²⟩ / (6 V kB t) with P = Σ qᵢ rᵢ over the ions; by
Nernst-Einstein, from the ions' self-diffusion; and their ratio, the ionicity.
The analysis is in :mod:`seamm_md_analysis.electric`.
"""

from collections import Counter
import importlib.resources
import logging
from pathlib import Path
import sys
import traceback
import warnings

import numpy as np
from seamm_md_analysis import open_md
from seamm_md_analysis.electric import (
    cell_dipole,
    classify_molecules,
    collective_msd,
    conductivity_from_msd,
    nernst_einstein,
)
from seamm_md_analysis.transport import fit_slope, msd_by_species, yeh_hummer_correction
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    preceding_viscosity,
    md_inputs,
    read_files,
    request_sampling,
    run_own_md,
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


def hill_formula(symbols):
    """The formula of a molecule in Hill order."""
    counts = Counter(symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(e for e in counts if e not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


class IonicConductivity(seamm.Node):
    """
    The non-graphical part of a Ionic Conductivity step in a flowchart.

    Attributes
    ----------
    subflowchart : seamm.Flowchart
        The MD to run, when this step runs its own MD.

    parameters : IonicConductivityParameters
        The control parameters for Ionic Conductivity.
    """

    def __init__(
        self,
        flowchart=None,
        title="Ionic Conductivity",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        logger.debug(f"Creating Ionic Conductivity {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Ionic Conductivity", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Ionic Conductivity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = ionic_conductivity_step.metadata
        self.parameters = ionic_conductivity_step.IonicConductivityParameters()
        self._file_handler = None
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

        self.subflowchart._parser = self.flowchart._parser

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do."""
        if P is None:
            P = self.parameters.values_to_dict()

        what = (
            "the ionic conductivity by Einstein-Helfand and Nernst-Einstein, and "
            "the ionicity"
        )

        source = P.get("md source", run_own_md)
        if source == run_own_md:
            text = f"Calculate {what}, from {P['nruns']} runs of the MD below.\n\n"
        else:
            text = f"Calculate {what}, from "
            if source == read_files:
                text += f"the MD in {P['md files']}.\n\n"
            else:
                text += f"{source} earlier in the flowchart.\n\n"
            short = True

        if not short:
            self.subflowchart.root_directory = self.flowchart.root_directory
            node = self.subflowchart.get_node("1").next()
            while node is not None:
                text += __(node.description_text()).__str__()
                text += "\n"
                node = node.next()

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Ionic Conductivity step."""
        next_node = super().run(printer)
        self._reset()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        _, configuration = self.get_system_configuration()
        if (
            P["md source"] == run_own_md
            and configuration is not None
            and configuration.periodicity != 3
        ):
            raise RuntimeError(
                "The Ionic Conductivity step needs a periodic system: a "
                "conductivity is a property of the bulk."
            )

        self.references.cite(
            raw=self._bibliography["Hansen2013"],
            alias="Hansen2013",
            module="ionic_conductivity_step",
            level=1,
            note="The Einstein-Helfand conductivity and the Nernst-Einstein relation.",
        )

        handlers = self._quiet_handlers()
        try:
            if P["md source"] == run_own_md:
                self._run_own_md(P)
            else:
                self.paths = md_inputs(self, P, name="all/position")
                printer.job(
                    __(
                        f"Analyzing {len(self.paths)} MD run"
                        f"{'s' if len(self.paths) > 1 else ''}: "
                        + ", ".join(str(p) for p in self.paths),
                        indent=self.indent + 4 * " ",
                    )
                )
        finally:
            self._restore_handlers(handlers)

        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        self.analyze(P=P, configuration=configuration)
        return next_node

    def _quiet_handlers(self):
        """Send only job-level output to job.out and the screen while MD runs."""
        result = []
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ) or isinstance(handler, logging.StreamHandler):
                result.append((handler, handler.level))
                handler.setLevel(printing.JOB)
        return result

    @staticmethod
    def _restore_handlers(handlers):
        for handler, level in handlers:
            handler.setLevel(level)

    def _run_own_md(self, P):
        """Run the MD in the sub-flowchart several times, keeping each run's file."""
        changes, missing = request_sampling(self.subflowchart, self.md_requirements(P))
        for change in changes:
            printer.job(__(f"Note: {change}.", indent=self.indent + 4 * " "))
        if len(missing) > 0:
            raise RuntimeError(
                "The MD in the sub-flowchart cannot provide what is needed: "
                + " ".join(missing)
            )

        first_node = self.subflowchart.get_node("1").next()
        if first_node is None:
            raise RuntimeError(
                "The sub-flowchart of the Ionic Conductivity step is empty."
            )
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()
        self.subflowchart.executor = self.flowchart.executor

        nruns = P["nruns"]
        fmt_run = f"0{len(str(nruns))}"
        for run in range(1, nruns + 1):
            run_dir = Path(self.directory) / f"run_{run:{fmt_run}}"
            run_dir.mkdir(parents=True, exist_ok=True)

            if self._file_handler is not None:
                self._file_handler.close()
                job.removeHandler(self._file_handler)
            path = run_dir / "Run.out"
            path.unlink(missing_ok=True)
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setLevel(printing.NORMAL)
            formatter = logging.Formatter(fmt="{message:s}", style="{")
            self._file_handler.setFormatter(formatter)
            job.addHandler(self._file_handler)

            self.flowchart.reset_visited()
            self.set_subids((*self._id, f"run_{run:{fmt_run}}"))

            node = first_node
            try:
                while node is not None:
                    node = node.run()
            except DeprecationWarning as e:
                printer.normal("\nDeprecation warning: " + str(e))
                traceback.print_exc(file=sys.stderr)
                traceback.print_exc(file=sys.stdout)
            except Exception as e:
                printer.job(f"Caught exception in run {run}: {str(e)}")
                with open(run_dir / "stderr.out", "a") as fd:
                    traceback.print_exc(file=fd)
                if "continue" in P["errors"]:
                    continue
                elif "exit" in P["errors"]:
                    break
                else:
                    raise
            else:
                self.paths.append(self._read_run(run, run_dir))

        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None

    def _read_run(self, run, run_dir):
        """The md.h5md of one run of the sub-flowchart with the positions."""
        paths = []
        for path in sorted(run_dir.glob("**/md.h5md")):
            with open_md(path) as md:
                if md.has_particles("position"):
                    paths.append(path)
        if len(paths) == 0:
            raise RuntimeError(f"There are no positions for run {run} in {run_dir}.")
        return paths[-1]

    # ------------------------------------------------------------------
    # The analysis
    # ------------------------------------------------------------------
    def _read(self, path, configuration, ion_atoms, ion_charges, ions, masses):
        """The collective charge displacement and the centres of mass of the ions
        in each frame, the interval between frames (ps), and the mean volume (Å³)
        and temperature (K) of one run."""
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
        return np.array(P), np.array(com), dt, float(np.mean(V)), T

    def analyze(self, indent="", P=None, configuration=None, **kwargs):
        """Analyze, print, plot and store the ionic conductivity."""
        if configuration is None:
            raise RuntimeError("There is no structure to analyze.")
        molecules = configuration.find_molecules(as_indices=True)
        symbols = configuration.atoms.symbols

        with open_md(self.paths[0]) as md:
            charges = md.charges("all")
            masses = md.masses("all")
        if charges is None:
            raise RuntimeError(
                f"'{self.title}': the MD has no fixed charges, as with a "
                "machine-learned forcefield, so there is no charge to follow."
            )
        charges = np.asarray(charges, dtype=float)
        masses = np.asarray(masses, dtype=float)

        classes = classify_molecules(charges, molecules)
        if len(classes["ions"]) == 0:
            raise RuntimeError(f"'{self.title}': there are no ions.")
        ions = [np.array(molecules[m]) for m in classes["ions"]]
        names = [hill_formula([symbols[i] for i in m]) for m in ions]
        z = {name: int(classes["z"][m]) for name, m in zip(names, classes["ions"])}
        species = list(dict.fromkeys(names))
        groups = [[i for i, n in enumerate(names) if n == s] for s in species]
        counts = {s: len(g) for s, g in zip(species, groups)}
        ion_atoms = np.concatenate(ions)
        ion_charges = charges[ion_atoms]
        printer.normal(
            __(
                "The ions: "
                + ", ".join(f"{counts[s]} {s} ({z[s]:+d})" for s in species)
                + ".",
                indent=self.indent + 4 * " ",
            )
        )

        runs = [
            self._read(p, configuration, ion_atoms, ion_charges, ions, masses)
            for p in self.paths
        ]
        nruns = len(runs)
        start, end = P["fit start"], P["fit end"]
        T = float(np.mean([r[4] for r in runs]))
        V = float(np.mean([r[3] for r in runs]))
        dt = runs[0][2] or 1.0
        n_frames = sum(len(r[0]) for r in runs)

        # Several runs: average the curves and take the spread of the runs. One
        # run: blocks of it for the error.
        if nruns > 1:
            pieces = [(r[0], r[1]) for r in runs]
        else:
            n_blocks = P["number of blocks"] if P["md source"] != run_own_md else 1
            n_blocks = max(1, min(int(n_blocks), len(runs[0][0]) // 20))
            pieces = list(
                zip(
                    np.array_split(runs[0][0], n_blocks),
                    np.array_split(runs[0][1], n_blocks),
                )
            )
        n = min(len(p[0]) for p in pieces)
        ts = np.arange(n) * dt
        msds = [collective_msd(Pc[:n]) for Pc, _ in pieces]
        sigmas = [
            conductivity_from_msd(m, dt, V, T, start=start, end=end)["sigma"]
            for m in msds
        ]
        D_pieces = []
        species_msds = []
        for _, com in pieces:
            msd, _ = msd_by_species(com[:n], groups)
            species_msds.append(msd)
            D_pieces.append(
                [
                    fit_slope(m[:, 3], ts.tolist(), start=start, end=end).slope
                    / 6
                    * 1.0e-8
                    for m in msd
                ]
            )

        # The value from the whole run(s); the error from the pieces
        if nruns > 1:
            mean_msd = np.mean(msds, axis=0)
            fit = conductivity_from_msd(mean_msd, dt, V, T, start=start, end=end)
            plot_msd, plot_ts = mean_msd, ts
        else:
            full = collective_msd(runs[0][0])
            full_ts = np.arange(full.size) * dt
            fit = conductivity_from_msd(full, dt, V, T, start=start, end=end)
            plot_msd, plot_ts = full, full_ts
        sigma = fit["sigma"]
        stderr = (
            float(np.std(sigmas, ddof=1) / np.sqrt(len(sigmas)))
            if len(sigmas) > 1
            else None
        )
        D = dict(zip(species, np.mean(D_pieces, axis=0).tolist()))
        D_err = (
            dict(
                zip(
                    species,
                    (
                        np.std(D_pieces, axis=0, ddof=1) / np.sqrt(len(D_pieces))
                    ).tolist(),
                )
            )
            if len(D_pieces) > 1
            else None
        )

        # The Yeh-Hummer correction of the self-diffusion coefficients
        notes = []
        D_corrected = None
        if P["Yeh-Hummer correction"] != "no":
            record = preceding_viscosity(self)
            if record is None:
                notes.append(
                    "The Yeh-Hummer correction is skipped because no Viscosity step "
                    "ran before this one."
                )
            else:
                L = V ** (1 / 3) * 1.0e-10
                correction = yeh_hummer_correction(T, L, record["value"])
                D_corrected = {s: d + correction for s, d in D.items()}
                notes.append(
                    "The self-diffusion coefficients are corrected for the size of "
                    f"the cell (Yeh-Hummer, +{correction * 1e9:.3g}e-9 m^2/s) with the "
                    f"viscosity from step {record['id']} ({record['step']})."
                )
        sigma_ne = nernst_einstein(D_corrected or D, counts, z, V, T)
        ionicity = sigma / sigma_ne if sigma_ne > 0 else None

        table = {"Property": [], "Value": [], "±": [], "Units": []}

        def row(name, value, err, units, fmt=".4g"):
            table["Property"].append(name)
            table["Value"].append("" if value is None else f"{value:{fmt}}")
            table["±"].append("" if err is None else f"{err:.2g}")
            table["Units"].append(units)

        row("σ (Einstein-Helfand)", sigma, stderr, "S/m")
        row("σ (Nernst-Einstein)", sigma_ne, None, "S/m")
        row("Ionicity σ/σ_NE", ionicity, None, "", fmt=".3f")
        for s in species:
            row(
                f"D({s})",
                D[s] * 1e9,
                None if D_err is None else D_err[s] * 1e9,
                "1e-9 m^2/s",
            )
            if D_corrected is not None:
                row(f"D({s}), corrected", D_corrected[s] * 1e9, None, "1e-9 m^2/s")
        text = tabulate(
            table, headers="keys", tablefmt="simple_outline", disable_numparse=True
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
        how = (
            f"the standard error over the {nruns} runs"
            if nruns > 1
            else (
                f"from {len(pieces)} blocks of the run"
                if len(pieces) > 1
                else "no error estimate: use several runs or blocks"
            )
        )
        printer.normal(
            __(
                f"From {n_frames} frames of {nruns} run{'s' if nruns > 1 else ''}, "
                f"every {dt:.3g} ps; the MSDs are fitted from {start:.0%} to {end:.0%} "
                f"of their length; the error is {how}. The collective displacement "
                "of the charges is noisy: long runs, or several, are needed.",
                indent=8 * " ",
            )
        )
        for note in notes:
            printer.normal(__(note, indent=8 * " "))
        printer.normal("")

        self._plot(plot_ts, plot_msd, fit, ts, species, species_msds)

        self.model = self._find_model()
        data = {
            "sigma": sigma,
            "sigma_NE": sigma_ne,
            "D": D,
            "T": T,
            "V": V,
            "nruns": nruns,
            "n_frames": n_frames,
        }
        if stderr is not None:
            data["sigma,stderr"] = stderr
        if ionicity is not None:
            data["ionicity"] = ionicity
        if D_corrected is not None:
            data["D, corrected"] = D_corrected
        self.store_results(configuration=configuration, data=data)

    def _plot(self, ts, msd, fit, species_ts, species, species_msds):
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
        mean = np.mean([[m[:, 3] for m in run] for run in species_msds], axis=0)
        for name, m in zip(species, mean):
            plot2.add_trace(
                x_axis=x2,
                y_axis=y2,
                name=name,
                x=list(species_ts),
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
        periodic cell."""
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
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed."""
        source = self._static("md source")
        if source is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze."""
        if self._static("md source") == run_own_md:
            return
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: a conductivity is a property of the bulk."
            )

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
