#!/usr/bin/env python

import rdkit, rdkit.Chem as rdkit
import stk
from random import shuffle
import subprocess as sp
import os, errno, shutil, itertools
from utils import *
from joblib import Parallel, delayed

class Session:

    """
    Builds a session, in which polymer structures may be constructed,
    optimized, have their properties calculated (IP, EA & S0->S1 transiton).
    Lists of monomers may be provided, co-polymer combinations of which
    may be automatically screened.

    parameters
    ----------

    name : :class:`str`
        Name of the session. Essentially the name of the directory into which
        all results will be placed.

    n_repeat : :class:`int`
        Number of repeat units that will be used to build a polymer.

    n_confs : :class:`int`
        Number of conformers to embed within conformer search.

    monomers_file : :class:`str`, optional
        A file containing a list of monomer units to be screened. This is
        required for screening with Session.screen().

    solvent : :class:`str` (default = ``None``)
        Solvent to be applied within calculations. Solvent effects are
        applied using an implicit model.

    random_select : :class:`bool` (default = ``False``)
        Randomly select co-polymer combinations to screen from 'monomers_file'.
        This may be useful if one requires randomly sampled compositions from
        the overall co-polymer space.

    returns
    -------
    str : :class:`str`
        A description of the current session.

    """

    def __init__(self, name, monomers_file=None, solvent=None):
        self.session_name = name
        self.solvent_info = []

        os.makedirs(self.session_name)

        with open(monomers_file) as f:
            self.monomers = [line.split() for line in f]

        if solvent is not None:
            if solvent not in valid_solvents:
                raise Exception('Invalid solvent choice. Valid solvents:',
                [i for i in valid_solvents])
            else:
                self.solvent_info = ['-gbsa', solvent]


    def screen(self, n_repeat, seq, nconfs, random_select=False):

        with open(self.session_name+'/'+'screening-output', 'w') as output:
            output.write(output_header)

        compositions = self.get_polymer_compositons(random_select)
        results = Parallel(n_jobs=15)(delayed(self.screening_protocol)(
                        n_repeat, seq, nconfs, random_select, composition)
                        for composition in compositions)


    def screening_protocol(self, n_repeat, seq, nconfs, random_select, composition):
        smiles1, id1 = composition[0][1], composition[0][0]
        smiles2, id2 = composition[1][1], composition[1][0]
        name = '{}-{}'.format(id1, id2)

        try:
            os.makedirs(self.session_name+'/'+name)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        with cd(self.session_name+'/'+name):
            try:
                polymer = self.generate_polymer(smiles1, smiles2, n_repeat, seq, name)
                conf, E = self.conformer_search(polymer, nconfs)
                E_xtb = self.xtb_opt(polymer)
                vip, vea = self.xtb_calc_potentials(polymer)
                gap, f = self.stda_calc_excitation(polymer)
                property_log(id1, id2, smiles1, smiles2, vip, vea, gap, f)

            except Exception as e:
                error_log(id1, id2, smiles1, smiles2, e)


    def get_polymer_compositons(self, random_select):
        homopolymers = [[i, i] for i in self.monomers]
        copolymers = list(itertools.combinations(self.monomers, 2))
        self.compositions = homopolymers + copolymers
        self.compositions.sort(key=lambda x: int(x[0][0]))

        if random_select:
            shuffle(self.compositions)

        return self.compositions


    def generate_polymer(self, smiles1, smiles2, n_repeat, seq, name):
        # initialise and embed rdkit mol objects
        a = rdkit.AddHs(rdkit.MolFromSmiles(smiles1))
        rdkit.AllChem.EmbedMolecule(a, rdkit.AllChem.ETKDG())
        b = rdkit.AddHs(rdkit.MolFromSmiles(smiles2))
        rdkit.AllChem.EmbedMolecule(b, rdkit.AllChem.ETKDG())

        # initialise stk StructUnit2 objects
        A = stk.StructUnit2.rdkit_init(a, "bromine")
        B = stk.StructUnit2.rdkit_init(b, "bromine")

        # construct polymer
        polymer = stk.Polymer([A,B], stk.Linear(seq, [0,0], n=n_repeat), name=name)
        stk.rdkit_ETKDG(polymer)

        return polymer


    def conformer_search(self, polymer, nconfs):

        mol = polymer.mol
        name = polymer.name
        confs = rdkit.AllChem.EmbedMultipleConfs(mol, nconfs, rdkit.AllChem.ETKDG())
        rdkit.SanitizeMol(mol)

        lowest_energy = 10**10
        for conf in confs:
            ff = rdkit.AllChem.MMFFGetMoleculeForceField(mol, rdkit.AllChem.MMFFGetMoleculeProperties(mol), confId=conf)
            ff.Initialize()
            energy = ff.CalcEnergy()

            if energy < lowest_energy:
                lowest_energy = energy
                lowest_conf = conf

        rdkit.MolToMolFile(mol, name+'.mol', confId=lowest_conf)

        return lowest_conf, lowest_energy


    def xtb_opt(self, polymer):
        name = polymer.name
        molfile = '{}.mol'.format(name)
        xyzfile = '{}.xyz'.format(name)
        sp.call(['babel', molfile, xyzfile])

        # Optimise, extract total & solv. energy
        calc_params = ['xtb', xyzfile, '-opt'] + self.solvent_info
        output = run_calc(calc_params)
        shutil.copy('xtbopt.xyz', '{}-opt.xyz'.format(name))


    def xtb_calc_potentials(self, polymer):
        name = polymer.name
        xyzfile = '{}-opt.xyz'.format(name)

        # calculate and extract IP
        calc_params = ['xtb', xyzfile, '-vip'] + self.solvent_info
        output = run_calc(calc_params)
        vip = output[output.find('delta SCC IP'):].split()[4]

        # calculate and extract EA
        calc_params = ['xtb', xyzfile, '-vea'] + self.solvent_info
        output = run_calc(calc_params)
        vea = output[output.find('delta SCC EA'):].split()[4]

        return vip, vea


    def stda_calc_excitation(self, polymer):
        name = polymer.name
        xyzfile = '{}-opt.xyz'.format(name)

        # calculate xtb wavefunction
        calc_params = ['xtb', xyzfile]
        run_calc(calc_params)

        # calculate excitations, extract S0 -> S1, extract f
        calc_params = ['stda', '-xtb', '-e', '8']
        output = run_calc(calc_params)
        gap = output[output.find('excitation energies'):].split()[13]
        f = output[output.find('excitation energies'):].split()[15]

        return gap, f


    def __str__(self):
        string = 'Session name: ' + self.session_name + '\n'

        return string
