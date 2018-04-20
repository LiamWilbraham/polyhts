
POLYHTS
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
* rdkit     http://www.rdkit.org/\n
* stk       https://github.com/supramolecular-toolkit/stk
* openbabel http://openbabel.org/wiki/Main_Page
* xtb       https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/xtb/xtb
* stda      https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/stda/stda

Input parameters for high throughput screening procedure
--------------------------------------------------------

    nconfs           : [int] Number of confomers generated in stochastic confomer search
                       By default, this is set to 1000

    length           : [int] Number of monomers in oligomer model (default = 8)

    xtb_dir          : [str] Full path to directory containing the executable for the xtb program 

    num_cores        : [int] Number of cores used in (embarrassingly) parallel polymer screening

    candidate_list   : [str] Name of a file containing a list of copolymer candidates. The candidate list file
                       must be written in the following format:

                       Unit 1    Unit 2    ID
                       smiles1   smiles2   XXXX

                       where smiles1 and smiles2 are SMILES strings that represent the monomer
                       units that will make up a copolymer with ID number XXXX.     

                       IMPORTANT - monomers should contain bromine atoms where they are to be connected
                       to adjacent monomer units.       

    solvent          : [str] Choose the solvent model parameters to be used. Available parameters are within 
                       the gbsa solvation model in GFN-xTB (acetone, acetonitrile, benzene, chcl3, cs2,
                       dmso, ether, h2o, methanol, thf, toluene). If no solvent model is to be included,
                       simply set this parameter to 'none'
 
    intensity_cutoff : [float] Excited states with oscillator strengths below this value will be rejected. This 
                       is to exclude excited states that will effectively not absorb light and thus 
                       contribute to the observed optical gap

    Linear calibration parameters :
                      - [float] ip_intercept      
                      - [float] ip_slope          
                      - [float] ip_intercept      
                      - [float] ea_slope          
                      - [float] opt_gap_intercept 
                      - [float] opt_gap_slope     
        
                      By default these are set to 1 and 0 for slopes and intercepts, respectively (i.e. no calibration)


Output files and format
-----------------------

    Results are formatted in a file named 'screened-polymers.dat':

    ID    IP      EA      Gap     F        E Solv.   Unit1                                       Unit2          
    0001  5.7464  3.5208  3.3690  7.0874   -3.0087   C1=CC2=C(C=C1Br)S(=O)(=O)C3=C2C=CC(=C3)Br   C1=C(SC(=N1)Br)Br                                           

* ID = polymer ID number
* IP = ionisation potential (in V)
* EA = electron affinity (in V)
* GAP = optical gap (in eV)
* F = oscillator strength
* E. Solv = solvation free energy (in eV)
* Unit1 = smiles string representing the first monomer unit
* Unit2 = smiles string representing the second monomer unit
