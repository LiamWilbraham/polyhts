#!/usr/bin/env python

import rdkit, rdkit.Chem as rdkit
import stk
from random import shuffle
import subprocess as sp
import os, errno, shutil, itertools
from utils import *

class Session:

    """
    Builds a session, in which polymer structures may be constructed,
    optimized, have their properties calculated (IP, EA & S0->S1 transiton).
    Lists of monomers may be provided, co-polymer combinations of which
    may be automatically screened.
    """

    def __init__(self, name, n_repeat, nconfs, A_monomers_file, solvent=None, random_select=False):
        self.name = name
        self.n_repeat = n_repeat
        self.nconfs = nconfs
        self.random_select = random_select
        self.solvent_info = []

        with open(A_monomers_file) as f:
            self.A_monomers = [line.split() for line in f]

        if solvent is not None:
            if solvent not in valid_solvents:
                raise Exception('Invalid solvent choice. Valid solvents:',
                [i for i in valid_solvents])
            else:
                self.solvent_info = ['-gbsa', solvent]


    def screen(self):
        with open('screening-output', 'w') as output:
            output.write(output_header)

        compositions = self.get_polymer_compositons()
        for composition in compositions:
            smiles1, id1 = composition[0][1], composition[0][0]
            smiles2, id2 = composition[1][1], composition[1][0]
            name = '{}-{}'.format(id1, id2)

            try:
                os.makedirs(name)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

            with cd(name):
                try:
                    polymer = self.generate_polymer(smiles1, smiles2)
                    conf, E = self.conformer_search(polymer, name)
                    E_xtb = self.xtb_opt(name)
                    vip, vea = self.xtb_calc_potentials(name)
                    gap, f = self.stda_calc_excitation(name)
                    property_log(id1, id2, smiles1, smiles2, vip, vea, gap, f)

                except Exception as e:
                    error_log(id1, id2, smiles1, smiles2, e)


    def get_polymer_compositons(self):
        homopolymers = [[i, i] for i in self.A_monomers]
        copolymers = list(itertools.combinations(self.A_monomers, 2))
        self.compositions = homopolymers + copolymers
        self.compositions.sort(key=lambda x: int(x[0][0]))

        if self.random_select:
            random.shuffle(self.compositions)

        return self.compositions


    def generate_polymer(self, smiles1, smiles2):
        # initialise and embed rdkit mol objects
        a = rdkit.AddHs(rdkit.MolFromSmiles(smiles1))
        rdkit.AllChem.EmbedMolecule(a, rdkit.AllChem.ETKDG())
        b = rdkit.AddHs(rdkit.MolFromSmiles(smiles2))
        rdkit.AllChem.EmbedMolecule(b, rdkit.AllChem.ETKDG())

        # initialise stk StructUnit2 objects
        A = stk.StructUnit2.rdkit_init(a, "bromine")
        B = stk.StructUnit2.rdkit_init(b, "bromine")

        # construct polymer
        polymer = stk.Polymer([A,B], stk.Linear("AB", [0,0], n=self.n_repeat))
        stk.rdkit_ETKDG(polymer)

        return polymer


    def conformer_search(self, polymer, name):

        mol = polymer.mol
        confs = rdkit.AllChem.EmbedMultipleConfs(mol, self.nconfs, rdkit.AllChem.ETKDG())
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


    def xtb_opt(self, name):
        molfile = '{}.mol'.format(name)
        xyzfile = '{}.xyz'.format(name)
        sp.call(['babel', molfile, xyzfile])

        # Optimise, extract total & solv. energy
        calc_params = ['xtb', xyzfile, '-opt'] + self.solvent_info
        output = run_calc(calc_params)
        shutil.copy('xtbopt.xyz', '{}-opt.xyz'.format(name))


    def xtb_calc_potentials(self, name):
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


    def stda_calc_excitation(self, name):
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
