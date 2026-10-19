.. _user-guide:

**********
User Guide
**********

The Ionic Conductivity step calculates the conductivity of an electrolyte from
equilibrium molecular dynamics, two ways:

* **Einstein-Helfand**, from the collective displacement of the ions' charges,

  .. math::

     \sigma = \lim_{t \to \infty} \frac{\langle \lVert \mathbf{P}(t) -
              \mathbf{P}(0) \rVert^2 \rangle}{6 V k_B t},
     \qquad \mathbf{P} = \sum_{\text{ions}} q_i \mathbf{r}_i

  which includes the correlations between the ions (ion pairing lowers it);
* **Nernst-Einstein**, from the self-diffusion coefficient Dₛ of each ion,
  σ_NE = e² / (V k_B T) Σₛ Nₛ zₛ² Dₛ, which treats the ions as independent;

and their ratio σ/σ_NE, the **ionicity**.

The ions are the molecules with a non-zero integer charge, monatomic like Na⁺ or
molecular like BF₄⁻; each ion's self-diffusion is that of its centre of mass. The
charges are the forcefield's: a machine-learned forcefield has no fixed charges, so the
step stops with a message for such MD. If the charges of the ions were scaled in the MD
(LAMMPS's "Scale the charges of ions by", e.g. 0.8), the ions are still found from their
formal charges, and both σ and σ_NE use the charges the ions carried in the MD, zₛ
times the factor: they are the conductivities of the scaled-charge model, and the
ionicity compares like with like.

Where the MD comes from
=======================

The step analyzes MD that ran before it in the flowchart, or the ``md.h5md`` files of
earlier runs. For independent runs, put the MD in a Loop -- in parallel, with each
run starting from the same structure ("Discard the structure" in the MD sub-steps)
and new velocities -- and analyze "the last MD run of each loop iteration", the
default. A typical flowchart:

.. code-block:: text

    Model Chemistry (FF@...), Packmol, minimization
    NPT equilibration; Set Cell to its mean volume; NVT equilibration
    Loop over the replicas (parallel)
        LAMMPS: velocities, NVT decorrelation, NVT production
    Viscosity (optional, for the Yeh-Hummer correction)
    Ionic Conductivity

The MD must sample the atoms' positions ("Sample the atomic positions" in the LAMMPS
NVE, NVT and NPT sub-steps) at least every picosecond ("Longest interval between
frames"); they are unwrapped, so the displacements are continuous. Checking the
flowchart before it runs reports what is missing, with a fix that turns the sampling
on, and also warns when the analysis is inside the loop holding its MD, when the loop
changes the temperature or pressure (the runs are then not replicas of one state), or
when the thermostat disturbs the dynamics (Langevin, Andersen, or NVE): use
Nosé-Hoover or Bussi for the production runs.

Before the MD runs, the step warns if the system is small for the analysis (fewer
than about 200 molecules) or has fewer than about 10 cation-anion pairs: 1 M salt in
1000 atoms is only about 6 pairs, too few for useful statistics. The ions are found
from the forcefield's charges once the atoms are typed, or else from the formal
charges; the same warning is given with the results. A system without ions is not an
error: the step says that there is no ionic conductivity to calculate and the job
goes on.

The collective displacement of the charges converges slowly: each run should be
several nanoseconds, and several runs are better than one.

Fits and errors
===============

The mean square displacements are fitted to straight lines between "Fit from" and "Fit
to", fractions of their length (by default 10 % to 50 %): the early part is not yet
diffusive, and the late part has few time origins. The collective displacement is a
single, noisy quantity per run -- this is why conductivities from MD need long runs,
or several.

The errors are from independent runs, or from blocks of a single run ("Blocks for
errors from one run"; the value is then from the whole run). The table gives the 95%
confidence interval -- Student's t times the standard error -- and says so; the
stored results carry the standard error.

For an isotropic run -- a cubic cell with nothing directional in it -- σ and each
ion's D are single values, the averages of their x, y and z components within each
run, and a line for each checks that the components agree within their errors.
Otherwise the components are given too. The values of each run or block are in the
detailed output. The graph ``ionic_conductivity.graph`` shows the collective MSD with
its fit and the MSD of each ion.

**Finite size.** The self-diffusion coefficients, and so σ_NE, depend on the size of
the cell. With "Yeh-Hummer correction" the step corrects them with the viscosity from
a Viscosity step earlier in the flowchart, as the Diffusivity step does. The table
then gives D and σ_NE both in the box and corrected, and a note gives the correction,
the viscosity used (with its standard error) and the step it came from; the ionicity
uses the corrected σ_NE. Without a Viscosity step, the correction is skipped with a
note.

Results
=======

* σ (S/m) by Einstein-Helfand, with its error;
* σ_NE by Nernst-Einstein, in the box and, with the correction, corrected, and the
  ionicity σ/σ_NE, with their errors;
* the self-diffusion coefficient of each ion, in the box and corrected;
* the state of the runs: temperature, pressure, volume and density with their errors,
  and for runs at constant volume the equilibrium density at the target pressure when
  an NPT run before them gives the bulk modulus;
* the numbers of runs and frames.

σ, σ_NE and the ionicity can be stored as properties of the configuration, and all
the results in variables and tables.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
