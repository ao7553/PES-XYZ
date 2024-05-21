# PES-XYZ

Fortran code ``xyz_pes.f90`` will construct the molecular potential energy surface of a three-atom XYZ molecule. 
Program will read in a list of potential energy surface expansion parameters followed by a user-defined grid of 
internal coordinates in the order:

<pre>r(X-Y) (in Angstrom)   r(Y-Z) (in Angstrom)   alpha(X-Y-Z) (in degrees)</pre>

where ``r(X-Y)`` is the bond length between atoms X and Y, ``r(Y-Z)`` is the bond length between atoms Y and Z, and ``alpha(X-Y-Z)``
is the interbond angle going between <X-Y-Z.

The program will print the nuclear geometry and corresponding potential energy (in wavenumbers) in the output file.

The folder ``inputs/`` contains three potential energy surfaces:

``ocs_pes.inp`` for the molecule carbonyl sulphide (OCS). See the OCS molecular line list publication <https://doi.org/10.1093/mnras/stae1110> for further details.

``koh_pes.inp`` for the molecule potassium hydroxide (KOH). See the KOH molecular line list publication <https://doi.org/10.1093/mnras/staa4041> for further details.

``naoh_pes.inp`` for the molecule sodium hydroxide (NaOH). See the NaOH molecular line list publication <https://doi.org/10.1093/mnras/staa4041> for further details.


