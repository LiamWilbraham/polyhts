from contextlib import contextmanager
import os
import subprocess as sp
from rdkit import Chem

valid_solvents = ['acetone', 'acetonitrile', 'benzene', 'chcl3', 'cs2',
                  'dmso', 'ether', 'h2o', 'methanol', 'thf', 'toluene']

output_header = 'ID IP EA S0->S1 F ESolv Smiles\n'

def print_formatted_properties(name, vip, vea, gap, f, E_solv):
    str = '{} :'.format(name)
    str += 'IP (eV)= {}  '.format(vip)
    str += 'EA (eV)= {}  '.format(vea)
    str += 'S0->S1 (eV)= {}  '.format(gap)
    str += 'F (a.u.)= {}  '.format(f)
    if E_solv is not None:
        str += 'E Solv. (eV)= {}  '.format(E_solv)
    print(str)


def factorial(x):
    return (1 if x==0 else x * factorial(x-1))


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def run_calc(calc_params):
    p = sp.Popen(calc_params, stdout=sp.PIPE, encoding='utf8')
    output, _ = p.communicate()
    return output


def property_log(repeat, vip, vea, gap, f, E_solv):
    smiles_repeat = Chem.MolToSmiles(Chem.RemoveHs(repeat.mol), canonical=True)
    with open('../screening-output', 'a+') as output:
        output.write('{} {} {} {} {} {} {}\n'.format(
        repeat.name, vip, vea, gap, f, E_solv, smiles_repeat))


def error_log(permutation, monomers_dict, e):
    with open('../screening-output', 'a+') as output:
        string = ''
        for id in permutation:
            string += monomers_dict[id] + ' '
        string += str(e)
        output.write('ERROR: {}\n'.format(string))


def remove_junk():
    os.system('rm charges charges3 energy tda.dat wfn.xtb wbo xtbopt* xtbrestart')
