%nproc=4
%mem=5760MB
%chk=gaussian_7_scn.chk
#P MP2/aug-cc-pVTZ scf(maxcycle=200) opt geom(AddGIC)

Gaussian input

-1 1
N      0.000000000000000      0.000000000000000     -1.820943000000000
C      0.000000000000000      0.000000000000000     -0.628586000000000
S      0.805131426886938      0.000000000000000      0.824196875196085

A(1,2,3) Freeze
