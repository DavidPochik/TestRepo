# athena-accretion (mhd-winds-new-EOS)
My version of mhd-winds-new-EOS (https://github.com/mattjraives/mhd-winds), which is fairly old (dates back to mid 2021).
The only updated portions since then are the qw_EOS.cpp file (Newest version provided by Tejas Prasanna on March 29th, 2022) and my problem generator file (accretion.cpp and similar copies).

For the accretion problem, I specify inflow at the inner and outer boundaries.
Values for v and rho at the outer boundary are calculated by using v_out = -sqrt(GM/(2R_out)) and rho_out=Mdot/(4 pi r^2 |v_out|).
For these calculations, Solar mass is assumed to be 2e33 g. Maybe I need a more specific value here, since this might not agree with how Solar mass is defined in Athena++.
Density is fixed at the inner boundary, where the value for the inner density is borrowed from the final output of the previous run, i.e., the input text file.
Values for the pressure at the inner and outer boundaries are borrowed from the final output of the previous run too.
Pressure is held constant because temperature needs to be held constant at the boundaries when interpolated from pressure.
Otherwise, unphysical rises in pressure and temperature occur at the boundaries.

The .athdf output files from the previous run are provided under 'PreviousRunFinalOutputData_DIR'.

This repository was created to diagnose issues with setting the mass accretion rate at the outer boundary for the accretion problem.
For instance, if a mass accretion rate of Mdot=1.2 Msun/s is specified in athinput.accretion at the outer boundary by holding the mass density and velocity constant, the resulting value for Mdot throughout the simulation will be close to 1.2 Msun/s, but not exact.

The default setup for the athinput.accretion file is to perturb the neutrino luminosity (from 16e51 erg/s to 20e51 erg/s for both electron neutrinos and antineutrinos) of a steady-state 1D accretion solution that has Mdot=1.2 Msun/s and Lnu=16e51 erg/s while keeping neutrino energies fixed.

# athena-accretion (accretion_pscalars)
I've borrowed a version of athena-1 (pulled around July 18th, 2022) and adapted my accretion problem to its passive scalar features. I've made changes to the following files from athena-1 to get the code to compile/run:
	1) eos.hpp
        2) qw_eos.cpp
        3) mesh.cpp
        4) interp_table.cpp

Whenever I run the code, there is a small chance that the EoS fails (obtains negative quantities), but this doesn't happen everytime. In other words, the code will sometimes fail at run time even if I make no changes to the code itself. This is more evident when running with mpi.

Another (possibly related) issue that happens is that the time increment dt randomely becomes very small when I run the code. It should be of order 10^-6 s with my current setup, but sometimes it gets set to something much smaller, and this seems to happen at random everytime I run the code. If I attempt to run enough times, I'll end up getting the correct dt value. Note that the solutions I obtain from the succesful runs seem physically reasonable.


