#!/home/liam/anaconda3/bin/python3.6

import stk
import os
from joblib import Parallel, delayed
import rdkit, rdkit.Chem as rdkit

'''
--------
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
 [1] J. Chem. Inf. Model. 2015, 121, 2562-2574 
 [2] J. Chem. Theory Comput. 2017, 13, 1989-2009
 [3] Comput. Theor. Chem. 2014, 1040, 45-53 
 [4] J. Chem. Phys. 2016, 145, 054103

--------------------------------------------------------
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

'''
# Input Parameters 
#---------------------------------------------
nconfs = 500
length = 8
xtb_dir = '/home/liam/software/XTB'
num_cores = 20
candidate_list = 'candidate-list.dat' 
solvent = 'h2o'
intensity_cutoff = 0.5
ip_intercept = 0.
ip_slope = 1.
ea_intercept = 0.
ea_slope = 1.
opt_gap_intercept = 0.
opt_gap_slope = 1.
#---------------------------------------------

# Generates a polymer with MTK
def generate_polymer(smiles1, smiles2, length, sequence, name):
    
    a = rdkit.MolFromSmiles(smiles1)
    a = rdkit.AddHs(a)
    rdkit.AllChem.EmbedMolecule(a, rdkit.AllChem.ETKDG())

    b = rdkit.MolFromSmiles(smiles2)
    b = rdkit.AddHs(b)
    rdkit.AllChem.EmbedMolecule(b, rdkit.AllChem.ETKDG())

    A = stk.StructUnit2.rdkit_init(a, "bromine")
    B = stk.StructUnit2.rdkit_init(b, "bromine")

    polymer = stk.Polymer([A,B], stk.Linear(sequence, [1,1], n=int(length/2)), name=name)
    polymer.write(name+'.mol')

    return polymer

# Performs confomer search with ETKDG and ranks with MMFF94
def confomer_search(nconfs, name):
    
    polymer = rdkit.MolFromMolFile(name+'.mol',sanitize=False)
    polymer_smiles = rdkit.MolToSmiles(polymer)
    polymer_new = rdkit.MolFromSmiles(polymer_smiles)
    polymer_new = rdkit.AddHs(polymer_new)
    rdkit.AllChem.EmbedMolecule(polymer_new,rdkit.AllChem.ETKDG())
    rdkit.MolToMolFile(polymer_new,'polymer_new.mol',forceV3000=True)
    rdkit.AllChem.EmbedMolecule(polymer_new, rdkit.AllChem.ETKDG())
    rdkit.MolToMolFile(polymer_new, name+'.mol', forceV3000=True)

    ids = rdkit.AllChem.EmbedMultipleConfs(polymer_new, nconfs, rdkit.AllChem.ETKDG())

    energies = []
    confomers = []   
    for id in ids:
        conf_name = "{}-{}.mol".format(name.replace('.mol', ''), id)
        rdkit.MolToMolFile(polymer_new, conf_name, confId=id)

        ff = rdkit.AllChem.MMFFGetMoleculeForceField(polymer_new, rdkit.AllChem.MMFFGetMoleculeProperties(polymer_new), confId=id)
        ff.Initialize()
        energies.append(ff.CalcEnergy())
        confomers.append(name+'-'+str(id))

    compiled = sorted(zip(confomers, energies), key=lambda x: x[1])
    id = compiled[0][0]
    lowest = compiled[0][1]

    with open(name+'-rankings','w') as rankings:
        for x, y in compiled:
            rankings.write(x+'    '+str(y)+'\n')

    return lowest, id


# Optimisation of lowest energy confomer with GFN-xTB
def confomer_opt(lowest, id, name, solvent):

    os.system('babel '+id+'.mol '+id+'.xyz')

    if solvent == 'none':
        os.system('xtb '+id+'.xyz -opt > '+id+'-opt.out')

    else:
        os.system('xtb '+id+'.xyz -opt -gbsa '+solvent+' > '+id+'-opt.out')    

    os.system('mv xtbopt.xyz '+id+'-opt.xyz')

    return id


# Calculation of IP with GFN-xTB using IP/EA fitted parameters
def confomer_ip(id, name, solvent):

    if solvent == 'none':
        os.system('xtb '+id+'-opt.xyz -vip > '+id+'-ip.out')

    else:
        os.system('xtb '+id+'-opt.xyz -vip -gbsa '+solvent+' > '+id+'-ip.out')

    with open(id+'-ip.out','r') as ip_out:
        inlines = [line.strip() for line in ip_out]
        for line in inlines:
            if "delta SCC IP (eV)" in line:
                ip_xtb = float(line.split()[-1])

    return ip_xtb


# Calculation of EA with GFN-xTB using IP/EA fitted parameters
def confomer_ea(id, name, solvent):
   
    if solvent == 'none':
        os.system('xtb '+id+'-opt.xyz -vea > '+id+'-ea.out')

    else:
        os.system('xtb '+id+'-opt.xyz -vea -gbsa '+solvent+' > '+id+'-ea.out')

    with open(id+'-ea.out','r') as ea_out:
        inlines = [line.strip() for line in ea_out]
        for line in inlines:
            if "delta SCC EA (eV)" in line:
                ea_xtb = float(line.split()[-1])

    return ea_xtb


# Calculate optical gap and oscillator strength with sTDA-xTB
def stda(id, name, solvent, intensity_cutoff):

    if solvent == 'none':
        os.system('xtb '+id+'-opt.xyz -gbsa '+solvent+' > '+id+'-wfn.out')

    else:
        os.system('xtb '+id+'-opt.xyz > '+id+'-wfn.out')

    os.system(' cp wfn.xtb '+id+'-stda.wfn')
    os.system('stda -xtb -e 10 > '+id+'-stda.out')

    with open(id+'-stda.out','r') as stda_out:
        inlines = [line.strip() for line in stda_out]

        for line_number, line in enumerate(inlines):
            if "state    eV      nm       fL        Rv(corr)" in line:
                begin = line_number

        for i in range(1,10):
            s1_data = inlines[begin + i]
            opt_gap = float(s1_data.split()[1])
            osc_strength = s1_data.split()[3]
            if float(osc_strength) > intensity_cutoff:
                break

    return opt_gap, osc_strength


# Calculate pseudo-B3LYP IP values using linear regression model
def ip_linearfit_convert(ip_xtb, ip_intercept, ip_slope):

    ip_b3lyp = str((ip_xtb - ip_intercept) / ip_slope)

    return ip_b3lyp


# Calculate pseudo-B3LYP EA values using linear regression model
def ea_linearfit_convert(ea_xtb, ea_intercept, ea_slope):

    ea_b3lyp = str((ea_xtb - ea_intercept) / ea_slope)

    return ea_b3lyp


# Calculate pseudo-TD-B3LYP optical gap using linear regression model
def opt_gap_linearfit_convert(opt_gap, opt_gap_intercept, opt_gap_slope):

    opt_gap_b3lyp = str((opt_gap - opt_gap_intercept) / opt_gap_slope)

    return opt_gap_b3lyp


# Pull (water) solvation energy from optimised GFN-xTB structure
def solv_energy(id, name):

    with open(id+'-opt.out','r') as solv_out:
        inlines = [line.strip() for line in solv_out]
        for line in inlines:
            if "gsolv         :" in line:
                solvenergy = 27.2114*float(line.split()[2])

    return solvenergy


# Read in list of polymer candidates from candidate list
def read_candidates(candidate_list):

    try:
        with open(candidate_list, 'r') as candidates:
            candidate_list = [line.split() for line in candidates]        

    except FileNotFoundError:
        print('ERROR : no list of candidates provided')
        exit()

    return candidate_list


# Confomer screening process, performed in parallel 
def main(smiles1, smiles2, tag, nconfs, length, solvent, intensity_cutoff):

    name = 'Polymer-'+tag
    os.system('mkdir '+name)
    os.chdir(name)

    try:
        generate_polymer(smiles1, smiles2, length, "AB", name)
 
        lowest, id = confomer_search(nconfs, name)

        id = confomer_opt(lowest, id, name,solvent)
        ip = ip_linearfit_convert(confomer_ip(id, name, solvent), ip_intercept, ip_slope)[:7]
        ea = ea_linearfit_convert(confomer_ea(id, name, solvent), ea_intercept, ea_slope)[:7]
        opt_gap, osc_strength = stda(id, name, solvent, intensity_cutoff)
        opt_gap = opt_gap_linearfit_convert(opt_gap, opt_gap_intercept, opt_gap_slope)[:7]
    
        if solvent != 'none':
            e_solv = str(solv_energy(id, name))[:7]

            with open('../screened-polymers.dat', 'a+') as screened:
                screened.write('{0}{1}{2}{3}{4}{5}{6}{7}\n'.format(
                               tag.ljust(6), ip.ljust(8), ea.ljust(8), opt_gap.ljust(8), osc_strength.ljust(9), 
                               e_solv.ljust(10), smiles1.ljust(60), smiles2.ljust(60)))

        else:
            with open('../screened-polymers.dat', 'a+') as screened:
                screened.write('{0}{1}{2}{3}{4}{5}{6}{7}\n'.format(
                               tag.ljust(6), ip.ljust(8), ea.ljust(8), opt_gap.ljust(8), osc_strength.ljust(9),
                               ('-----').ljust(10), smiles1.ljust(60), smiles2.ljust(60)))


    except Exception as error:
        with open('../screened-polymers.dat', 'a+') as screened:
            screened.write('{0} ERROR: {1}\n'.format(tag.ljust(6), str(error).ljust(42)))

    os.chdir('../')


#--------------------------------------------------------------------------------------------------
# Run polymer screening in parallel with number of parallel threads = num_cores
if __name__ == "__main__":

    valid_solvent_list = ['acetone', 'acetonitrile', 'benzene', 'chcl3', 'cs2',
                          'dmso', 'ether', 'h2o', 'methanol', 'thf', 'toluene', 'none']

    error = True
    for solvents in valid_solvent_list:
        if solvent in solvents:
            error = False

    if error == True:
        print('ERROR: invalid solvent')
        exit()

    candidates = read_candidates(candidate_list)

    os.environ['XTBHOME']= xtb_dir

    with open('screened-polymers.dat', 'w') as screened:
        screened.write('ID    IP      EA      Gap     F        E Solv.   Unit1                                                       Unit2          \n')

    try:
        results = Parallel(n_jobs=num_cores)(delayed(main)(smiles1, smiles2, tag, nconfs, length, solvent, intensity_cutoff) for smiles1, smiles2, tag in candidates)

    except (ValueError):
        print('ERROR : perhaps the candidate list file is in an invalid format')
        exit()
#--------------------------------------------------------------------------------------------------
