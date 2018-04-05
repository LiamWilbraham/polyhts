#!/home/liam/anaconda3/bin/python3.6

import mtk
import os
from joblib import Parallel, delayed
import rdkit, rdkit.Chem as rdkit

# Generates a polymer with MTK
def generate_polymer(unit1, unit2, length, sequence, name, monomer_dir):

    A = mtk.StructUnit2(monomer_dir+unit1, "bromine")
    B = mtk.Unit2(monomer_dir+unit2, "bromine")
    polymer = mtk.Polymer([A,B], mtk.Linear(sequence, [1,1], n=int(length/len(sequence))), name=name)
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
def confomer_opt(lowest, id, name):

    os.system('babel '+id+'.mol '+id+'.xyz')
    os.system('xtb '+id+'.xyz -opt -gbsa h2o > '+id+'-opt.out')    
    os.system('mv xtbopt.xyz '+id+'-opt.xyz')

    return id


# Calculation of IP with GFN-xTB using IP/EA fitted parameters
def confomer_ip(id, name):

    os.system('xtb '+id+'-opt.xyz -vip -gbsa h2o > '+id+'-ip.out')

    with open(id+'-ip.out','r') as ip_out:
        inlines = [line.strip() for line in ip_out]
        for line in inlines:
            if "delta SCC IP (eV)" in line:
                ip_xtb = float(line.split()[-1])

    return ip_xtb


# Calculation of EA with GFN-xTB using IP/EA fitted parameters
def confomer_ea(id, name):
   
    os.system('xtb '+id+'-opt.xyz -vea -gbsa h2o > '+id+'-ea.out')

    with open(id+'-ea.out','r') as ea_out:
        inlines = [line.strip() for line in ea_out]
        for line in inlines:
            if "delta SCC EA (eV)" in line:
                ea_xtb = float(line.split()[-1])

    return ea_xtb


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


# Calculate optical gap and oscillator strength with sTDA-xTB
def stda(id, name):

    os.system('xtb '+id+'-opt.xyz -gbsa h2o > '+id+'-wfn.out')
    os.system(' cp wfn.xtb '+id+'-stda.wfn')
    os.system('stda -xtb -e 10 > '+id+'-stda.out')

    with open(id+'-stda.out','r') as stda_out:
        inlines = [line.strip() for line in stda_out]

        for line_number, line in enumerate(inlines):
            if "state    eV      nm       fL        Rv(corr)" in line:
                begin = line_number + 1

        s1_data = inlines[begin]
        opt_gap = float(s1_data.split()[1])
        osc_strength = s1_data.split()[3]

    return opt_gap, osc_strength


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

    with open(candidate_list, 'r') as candidates:
        candidate_list = [line.split() for line in candidates]        
        candidate_list = candidate_list[1:]

    return candidate_list


# Confomer screening process, performed in parallel 
def main(unit1, unit2, tag, nconfs, monomer_dir):

    name = 'Polymer-'+tag
    os.system('mkdir '+name)
    os.chdir(name)

    try:
        unit1 = 'monomer_'+unit1+'.mol'
        unit2 = 'monomer_'+unit2+'.mol'
        generate_polymer(unit1, unit2, 8, "AB", name, monomer_dir)
 
        lowest, id = confomer_search(nconfs, name)

        id = confomer_opt(lowest, id, name)
        ip = ip_linearfit_convert(confomer_ip(id, name), ip_intercept, ip_slope)[:7]
        ea = ea_linearfit_convert(confomer_ea(id, name), ea_intercept, ea_slope)[:7]
        opt_gap, osc_strength = stda(id, name)
        opt_gap = opt_gap_linearfit_convert(opt_gap, opt_gap_intercept, opt_gap_slope)[:7]
        e_solv = str(solv_energy(id, name))[:7]

        with open('../screened-polymers.dat', 'a+') as screened:
            screened.write('{0}  {1}  {2}  {3}  {4}  {5}  {6}{7}\n'.format(
                           unit1[-8:-4].ljust(10), unit2[-8:-4].ljust(10), tag.ljust(10), 
                           ip.ljust(10), ea.ljust(10), opt_gap.ljust(18), 
                           osc_strength.ljust(18), e_solv.ljust(5)))

    except (OSError, TypeError,NameError, ValueError, AttributeError):
        with open('../screened-polymers.dat', 'a+') as screened:
            screened.write('FAILURE '+tag+'\n')

    os.chdir('../')

#################################################################################
'''
Input parameters for high throughput screening procedure
----------------------------------------------------------

    nconfs      : Number of confomers generated in stochastic confomer search
                  By default, this is set to 1000

    monomer_dir : Directory containing .mol files of monomer units used to construct copolymers
                  By default, this is set to a directory 'monomers' within the working directory

    num_cores   : Number of cores used in (embarrassingly) parallel polymer screening

    Linear calibration parameters :
        - ip_intercept      
        - ip_slope          
        - ea_intercept      
        - ea_slope          
        - opt_gap_intercept 
        - opt_gap_slope     
        
        By default these are set to 1 and 0 for slopes and intercepts, respectively (i.e. no calibration)

'''

nconfs = 1000
monomer_dir = os.getcwd()+'/monomers/'
num_cores = 20
ip_intercept = 0. 
ip_slope = 1.
ea_intercept = 0.
ea_slope = 1.
opt_gap_intercept = 0.
opt_gap_slope = 1.

# Run polymer screening in parallel with number of parallel threads = num_cores
if __name__ == "__main__":
    candidates = read_candidates('candidate_list.dat')
    with open('screened-polymers.dat', 'w') as screened:
        screened.write('Unit1       Unit2       ID          IP (eV)     EA (eV)     Opt. Gap (eV)       Osc. Strength     Solvation Energy(eV)\n')
    results = Parallel(n_jobs=num_cores)(delayed(main)(unit1, unit2, tag, nconfs, monomer_dir) for unit1, unit2, tag in candidates)

#################################################################################
