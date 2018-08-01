
POLYHTS :sunny:
--------
A high throughput screening script for computing optoelectronic and redox properties of organic polymers.
Starting with the structures of candidate monomer units, arbitrary co-polymers can be constructed and 
screened for light absorption, ionisation potentials, electron affinities, and solvation free energies. 

Before the calculation of the above properties, a confomer search is conducted stochastically using the
Experimental Torsion Distance Geometry with additional basic Knowledge (ETKDG) [1] approach. On the resulting
lowest-energy confomer, structure optimisation and electronic properties are calculated using a family of 
semi-empirical methods based on the GFN-xTB density functional tight binding method proposed by Grimme [2-4].
By default, copolymers are modelled using a linear chain of 12 monomer units. Solvation effects can be included
using the gbsa solvation model associated with GFN-xTB.

Additonally, one may provide parameters to calibrate the computed data. For example, for a particular problem
you might want to calibrate the semi-empirically derived properties to those computed using a higher-level 
electronic structure method. This can be done automatically by providing linear regression fit parametes
(i.e. slopes and intercepts) to the script (see below).

references
----------
* [1] J. Chem. Inf. Model. 2015, 121, 2562-2574 
* [2] J. Chem. Theory Comput. 2017, 13, 1989-2009
* [3] Comput. Theor. Chem. 2014, 1040, 45-53 
* [4] J. Chem. Phys. 2016, 145, 054103

requirements
------------
* rdkit     http://www.rdkit.org/
* stk       https://github.com/supramolecular-toolkit/stk
* openbabel http://openbabel.org/wiki/Main_Page
* xtb       https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/xtb/xtb
* stda      https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/stda/stda
