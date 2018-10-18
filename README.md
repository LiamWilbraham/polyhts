
# `polyhts` :sunny:
Python library for ultra-fast structural optimisation and property calculations of organic co-polymers.

`polyhts` takes base functionality from:

* `rdkit` : transforms SMILES strings into monomer building blocks
* `stk` : constructs linear co-polymer structures from these building blocks
* `xtb` : optimises geometries, calculates properties

Combining the above, `polyhts` can be used for high-speed, accurate screening of
organic co-polymer compositions, both for property screening and quick, exploratory
calculations. Currently supported properties are:

* Ionisation potentials (IP)
* Electron affinities (EA)
* Excitation energies & oscillator strengths
* Solvation free energies

Though, in principle, any property that can be calculated via `xtb` or `stda`
can be obtained.

## Functionality
`polyhts` calculations start by defining a `Session`, in which information like
co-polymer repeat unit length, number of repeat units that will be used to
construct a polymer chain, number of conformers to be explored and solvent type
are specified.

For example, we can start a `Session` in which we will construct polymers with 4
repeat units (each of which are 4 monomers long) and explore 100 conformers while
applying an implicit solvent model for benzene:

```python
session = polyhts.Session('my_session', 2, 4,  100, solvent='benzene')  
```
### 1. Combinatorial Screening
Within this session, we can screen all combinations of pre-supplied monomer unit
SMILES from a text file:
```python
session.screen('smiles-list.txt', nprocs=20)      
```
where `smiles-list.txt` has the format:
```
0001 smiles1
0002 smiles2
0003 smiles3
 .     .
 .     .
 .     .
```
Note that, not only will all compositions of monomers be screened, but all
permutations of a given compositions as well (e.g. AABB as well as ABAB).

### 2. Fix one monomer, screen possible co-monomers (under construction)

### 3. Just one Polymer
We can also calculate properties for a single co-polymer, where we supply a list
of smiles explicitly, where the length of this list is equivalent to the
repeat unit length:

```python
session.calc_polymer_properties(['c1c(Br)cc(Br)cc1', 'c1c(Br)cc(Br)cc1'], 'polymer-name')  
```
Following the `stk` documentation, `Br` atoms are places within SMILES strings
where monomer units are to be connected to one another.

## Installation

## references
* [1] J. Chem. Inf. Model. 2015, 121, 2562-2574  
* [2] J. Chem. Theory Comput. 2017, 13, 1989-2009
* [3] Comput. Theor. Chem. 2014, 1040, 45-53
* [4] J. Chem. Phys. 2016, 145, 054103

## requirements
* rdkit     http://www.rdkit.org/
* stk       https://github.com/supramolecular-toolkit/stk
* openbabel http://openbabel.org/wiki/Main_Page
* xtb       https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/xtb/xtb
* stda      https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/stda/stda
