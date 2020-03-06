# Aurelia-force-and-torque
Script used to calculate forces and torques around Aurelia aurita jellies from surface pressures during turning maneuvers
(used in https://doi.org/10.1101/706762)


Reads in:

-images used for PIV processing

-velocity vector fields from PIV

-outlines of the jellyfish

-pressure field data from John Dabiri's queen2 (linked below)
http://dabirilab.com/software

Requires that the following Matlab scripts are loaded:

-interparc: https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc

-polygeom: https://www.mathworks.com/matlabcentral/fileexchange/319-polygeom-m

-cbrewer (optional): https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab

For background on force calculation methods, see:
https://doi.org/10.1371/journal.pone.0189225

For background on pressure calculation via queen2, see:
https://doi.org/10.1242/jeb.092767

