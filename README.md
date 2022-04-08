# athena-accretion
My version of mhd-winds-new-EOS (https://github.com/mattjraives/mhd-winds), which is fairly old (dates back to mid 2021).
The only updated portions since then are the qw_EOS.cpp file (Newest version provided by Tejas Prasanna on March 29th, 2022) and my problem generator file (accretion.cpp and similar copies).

For the accretion problem, I specify inflow at the inner and outer boundaries.
Values for v and rho at the outer boundary are calculated by using v_out = -sqrt(GM/(2R_out)) and rho_out=Mdot/(4 pi r^2 |v_out|).
For these calculations, Solar mass is assumed to be 2e33 g. Maybe I need a more specific value here.
Density is fixed at the inner boundary, where the value for the inner density is borrowed from the final output of the previous run, i.e., the input text file.
Values for the pressure at the inner and outer boundaries are borrowed from the final output of the previous run too.
Pressure is held constant because temperature needs to be held constant at the boundaries when interpolated from pressure.
Otherwise, unphysical rises in pressure and temperature occur at the boundaries.

The .athdf output files from the previous run are provided under 'PreviousRunFinalOutputData_DIR'.

This repository was created to diagnose issues with setting the mass accretion rate at the outer boundary for the accretion problem.
For instance, if a mass accretion rate of Mdot=1.2 Msun/s is specified in athinput.accretion at the outer boundary by holding the mass density and velocity constant, the resulting value for Mdot throughout the simulation will be close to 1.2 Msun/s, but not exact.

The default setup for the athinput.accretion file is to perturb the neutrino luminosity (from 16e51 erg/s to 20e51 erg/s for both electron neutrinos and antineutrinos) of a steady-state 1D accretion solution that has Mdot=1.2 Msun/s and Lnu=16e51 erg/s while keeping neutrino energies fixed.
