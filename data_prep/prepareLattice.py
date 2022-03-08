#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 12:18:42 2021

@author: hadfield

script to prepare the lattice for manual pharmacophore specification
"""

import sys
from rdkit import Chem
import numpy as np
import shutil
import os

sys.path.insert(0, '..')
from preprocessing import preprocessing




###Specify inputs#####

mol3DPath = sys.argv[1] #Location of the 3D conformer of the fragment
                        #Can also be a 3D conformer of a molecule of which the fragment is a substructure,
                        #in which case the fragment conformer will be extracted


fragSmilesInput = sys.argv[2] #location of the txt file containing the fragment smiles or a string containing the fragment SMILES
proteinLoc = sys.argv[3] #We're just going to copy the protein file into outStoreDir for easier use in the pymol script

outStoreDir = sys.argv[4] #Directory to store files in

#########################


if fragSmilesInput is not None:
    
    #Check whether args.fragment_smiles is a file
    if os.path.exists(os.path.expanduser(fragSmilesInput)):
        fragSmilesInput = os.path.abspath(os.path.expanduser(fragSmilesInput))
        smiles_file = True
    else:
        #Check that we can parse it with RDKit
        try:
            if Chem.MolFromSmiles(args.fragment_smiles).GetNumHeavyAtoms() > 0:
                smiles_file = False
            else: 
                raise ValueError('Fragment must have at least one heavy atom')
        except:
            
            raise ValueError("The supplied fragment_smiles doesn't appear to be a file and RDKit is unable to parse it as a molecule. Please check that you're providing a valid fragment")

    if smiles_file:
        with open(fragSmilesInput, 'r') as f:
            fragSmiles = f.read()
    
    else:
        fragSmiles = fragSmilesInput





prep = preprocessing()
mol3D = Chem.SDMolSupplier(mol3DPath)[0]


###Prepare fragment###
fc, evI, evp, fc2 = prep.preprocessFragment(fragSmiles, mol3D)

#Save fragment SDF
Chem.MolToMolFile(fc, f'{outStoreDir}/frag.sdf')

#Save constraint SDF (will need to be converted to Mol2 using obabel)
Chem.MolToMolFile(fc2, f'{outStoreDir}/constraint.sdf')

#Save fragment exit position
np.savetxt(f'{outStoreDir}/evp.txt', evp)


#Generate lattice
lattice = prep.getLattice(evp, size = 9, resolution = 1) #size controls the dimension of the lattice (in Angstroms), resolution controls how tightly packed the lattice points are (smaller = more tightly packed) 

#Create a series of atoms, positioned in the coordinates specified in the lattice
mols = []

for idx, coords in enumerate(lattice):
    mol = Chem.MolFromSmiles('I')
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, coords)

    conf.SetId(0)
    mol.AddConformer(conf)

    mols.append(mol)

#Write to file
w = Chem.SDWriter(f'{outStoreDir}/pharmacophoreLattice.sdf')
for m in mols:
    w.write(m)

#Copy protein into outStoreDir
shutil.copy(proteinLoc, f'{outStoreDir}/protein.pdb')
