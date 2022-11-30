# athena-accretion (helmholtz + passive scalars)
'athena_helmholtz_passive_scalars/' is a slightly modified version of the temperature_eos branch in athena_1 (pulled on 11/22/2022).
I've added a new helmholtz eos file called 'helmholtz_experimental', which is adapted to work with the accretion problem generator.

The files in 'athena_helmholtz_passive_scalars/' that differ from athena-1 include the following:

	1) configure.py

		- Changed 'helmholtz' -> 'helmholtz_experimental' in the conditional statement that enables a table EoS

	2) src/eos/eos.hpp

		- Added PresFromRhoT(...) and TFromRhoP(...) EoS functions to accommodate EoS calls in the accretion problem generator.

	3) src/eos/general/qw_eos.cpp

		- Added PresFromRhoT(...) and TFromRhoP(...) EoS functions with passive scalar capabilities.

		- Added a conditional statement to AsqFromRhoP(...) to prevent issues at runtime. This is a bit of a hack.

	4) src/utils/eos_table_class.cpp

		- Included helmholtz table options to accommodate the helm_table.dat filetype

