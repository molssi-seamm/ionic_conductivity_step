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
step stops with a message for such MD.

Where the MD comes from
=======================

As for the other analysis steps: the step runs its own MD a number of times, or
analyzes the MD earlier in the flowchart or the ``md.h5md`` files of earlier runs. The
MD must sample the atoms' positions at least every picosecond ("Longest interval
between frames"); they are unwrapped, so the displacements are continuous.

Fits and errors
===============

The mean square displacements are fitted to straight lines between "Fit from" and "Fit
to", fractions of their length (by default 10 % to 50 %): the early part is not yet
diffusive, and the late part has few time origins. The collective displacement is a
single, noisy quantity per run -- this is why conductivities from MD need long runs,
or several. The error is the standard error over the runs, or over blocks of a single
run. The graph ``ionic_conductivity.graph`` shows the collective MSD with its fit and
the MSD of each ion.

**Finite size.** The self-diffusion coefficients, and so σ_NE, depend on the size of
the cell. With "Yeh-Hummer correction" the step corrects them with the viscosity from
a Viscosity step earlier in the flowchart, as the Diffusivity step does, and says so;
without one, the correction is skipped with a note.

Results
=======

σ (S/m) and its standard error, σ_NE, the ionicity (which can be stored as properties
of the configuration), the self-diffusion coefficient of each ion (and corrected), the
mean temperature and volume, and the numbers of runs and frames.

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
