#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 11:46:08 2021

@author: hadfield
"""

from rdkit import Chem
from rdkit.Chem.Draw import MolToFile
import sys


def addDummyAtomToMol(mol, atomIdx):
    #Provide a molecule and atom index
    #returns a molecule with a dummy atom connected to that atom
    
    rwmol = Chem.RWMol(mol)
    dummyAtomIdx = rwmol.AddAtom(Chem.Atom(0)) #Add dummy atom and get its idx
    rwmol.AddBond(atomIdx, dummyAtomIdx, Chem.BondType.SINGLE)
    
    Chem.SanitizeMol(rwmol)
    return rwmol
    

def doCanonicalAtomRenumbering(mol):
    #Get canonical atom ordering and renumber 
    #return the renumbered molecule
    mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]#
    mol_renum = Chem.RenumberAtoms(mol, mol_neworder)
    
    return mol_renum


def addAtomIndices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i)


def getMolImageWithAtomNumbers(mol, imgFileName):
    
    mol = doCanonicalAtomRenumbering(mol)
    addAtomIndices(mol)
    MolToFile(mol, imgFileName, size = (500, 500))
    return 0


def getSMILESStringForSTRIFE(molSDF, imgFileName = 'fragmentWithAtomNumbering.png', smilesFileName = None):
    
    #Function which takes a molecule file as input and allows the user to select an exit vector. 
    #Outputs a SMILES string including the exit vector which can then be passed to STRIFE
    #Arguments:
    #molSDF - An SDF file containing the fragment
    #imgFileName - Location to save the png file containing the numbered fragment image
    #smilesFileName - Optionally specify a txt file to store the resulting SMILES string 
    #returns: SMILES string with a dummy atom denoting the exit vector.
    
    
    
    mol = Chem.MolFromMolFile(molSDF)
    
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol)) #Make mol2D so can visualise more easily
    mol = doCanonicalAtomRenumbering(mol)
    getMolImageWithAtomNumbers(mol, imgFileName)
    
    #Choose atom idx
    atom_idx = int(input(f'Please select an atom index for the fragment exit vector: An image showing the index associated with each atom can be found at {imgFileName}.'))
    molWithDummy = addDummyAtomToMol(mol, atom_idx)
    
    smiles = Chem.MolToSmiles(molWithDummy)
    smiles = reformatDummyAtomInSMILES(smiles)
    
    if smilesFileName is not None:
        with open(smilesFileName, 'w') as f:
            f.write(smiles)
    
    return smiles
    
def reformatDummyAtomInSMILES(smilesString):
    
    if '[*:' in smilesString:
        print('No dummy atom reformatting necessary. Returning input string')
        return smilesString
    else:
        smilesString = smilesString.replace('*', '[*:1]')
        return smilesString

if __name__=='__main__':
    

#Usage:
    #python specifyExitVector.py <fragmentSDF> <imageFileName> (<smilesFileName> (optionally))
    

    molSDF = sys.argv[1]
    imgFileName = sys.argv[2]
    
    if len(sys.argv) > 3:
        smilesFileName = sys.argv[3]
    else:
        smilesFileName = None
        
    smilesInputToSTRIFE = getSMILESStringForSTRIFE(molSDF, imgFileName, smilesFileName)
    
    
    
    print(f'You can use the following SMILES string as an input to STRIFE: {smilesInputToSTRIFE}')


    
    
    
    
    

