from contextlib import contextmanager
import os
import subprocess as sp

valid_solvents = ['acetone', 'acetonitrile', 'benzene', 'chcl3', 'cs2',
                  'dmso', 'ether', 'h2o', 'methanol', 'thf', 'toluene']

output_header = 'ID1\tID2\tIP\tEA\tS0->S1\tF\tESolv\tSmiles1\tSmiles2\n'

def print_formatted_properties(name, vip, vea, gap, f, E_solv):
    print(name, ':')
    print('IP (V) =', vip)
    print('EA (V) =', vea)
    print('S0 -> S1 (eV) =', gap)
    print('F (a.u.) =', f)
    if E_solv is not None:
        print('E Solv. =', E_solv)


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


def property_log(id1, id2, smiles1, smiles2, vip, vea, gap, f, E_solv):
    with open('../screening-output', 'a+') as output:
        output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        id1, id2, vip, vea, gap, f, E_solv, smiles1, smiles2))


def error_log(id1, id2, smiles1, smiles2, e):
    with open('../screening-output', 'a+') as output:
        output.write('{}\t{}\t{}\t{}\tERROR:{}\n'.format(id1, id2, smiles1, smiles2, e))


def remove_junk():
    os.system('rm charges charges3 energy tda.dat wfn.xtb wbo xtbopt* xtbrestart')
