=======
History
=======

2026.10.16 -- First release: the ionic conductivity
    * The ionic conductivity of electrolytes from MD with fixed charges, by
      Einstein-Helfand from the collective displacement of the ions' charges and by
      Nernst-Einstein from their self-diffusion, and the ionicity; errors from
      independent runs or blocks; the Yeh-Hummer correction of the self-diffusion
      with the viscosity from a Viscosity step.
    * Runs its own MD or analyzes the MD earlier in the flowchart, and checks before
      the job runs that the MD samples the positions often enough.
