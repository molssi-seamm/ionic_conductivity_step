=======
History
=======
2026.10.19 -- Scaled ion charges
    * The ions of MD with scaled charges (LAMMPS's 'Scale the charges of ions by', e.g.
      0.8) were not recognized, since their charges are not integers. They are now
      found from their formal charges, and the Nernst-Einstein conductivity uses the
      charges the ions carried in the MD, like the Einstein-Helfand one, so that the
      ionicity compares like with like.

2026.10.17.1 -- Analyze MD run before the step; labelled errors; no ions skipped
    * The step analyzes MD run before it in the flowchart -- the runs of a Loop around
      the MD, as independent replicas ("the last MD run of each loop iteration", the
      default) -- or md.h5md files; it no longer runs its own MD in a sub-flowchart.
      Older flowcharts that stored its settings for that can still be read; the check
      before the job runs says what replaces them.
    * The check before the job runs also warns when the analysis is inside the loop
      holding its MD, when the loop changes the temperature or pressure, and when the
      thermostat does not suit the analysis, with fixes; the step warns before the MD
      when the system is small.
    * Errors are the 95% confidence interval, Student's t for the number of runs or
      blocks, and the output says so; stored errors are standard errors.
    * The state of the runs -- temperature, pressure, volume and density with their
      errors, and the equilibrium density where an NPT run before them gives the bulk
      modulus -- is reported and stored. The values of each run or block are in
      step.out.
    * Bugfix: a system without ions stopped the job; it is now reported and skipped.
    * A warning before the MD runs when there are fewer than about 10 ion pairs.
    * σ and the ions' diffusion coefficients are averages of their x, y and z components,
      with a check that they agree; with the Yeh-Hummer correction, the values in the
      box and corrected are both given, with the viscosity used.


2026.10.17 -- First release: the ionic conductivity
    * The ionic conductivity of electrolytes from MD with fixed charges, by
      Einstein-Helfand from the collective displacement of the ions' charges and by
      Nernst-Einstein from their self-diffusion, and the ionicity; errors from
      independent runs or blocks; the Yeh-Hummer correction of the self-diffusion
      with the viscosity from a Viscosity step.
    * Runs its own MD or analyzes the MD earlier in the flowchart, and checks before
      the job runs that the MD samples the positions often enough.
