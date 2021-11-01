#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 11:46:08 2021

@author: hadfield
"""

from rdkit import Chem
from rdkit.Chem.Draw import MolToFile
import sys
import argparse


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

def removeAtomIndices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(0)



def getMolImageWithAtomNumbers(mol, imgFileName):
    
    #mol = doCanonicalAtomRenumbering(mol)
    addAtomIndices(mol)
    MolToFile(mol, imgFileName, size = (500, 500))
    return 0


def getSMILESStringForSTRIFE(molSDF, imgFileName = 'fragmentWithAtomNumbering.png', returnSmiles = False, smilesFileName = None):
    
    #Function which takes a molecule file as input and allows the user to select an exit vector. 
    #Outputs a SMILES string including the exit vector which can then be passed to STRIFE
    #Arguments:
    #molSDF - An SDF file containing the fragment
    #imgFileName - Location to save the png file containing the numbered fragment image
    #smilesFileName - Optionally specify a txt file to store the resulting SMILES string 
    #returns: SMILES string with a dummy atom denoting the exit vector.
    
    
    
    mol = Chem.MolFromMolFile(molSDF)
    
    mol.RemoveAllConformers() #Make mol2D so can visualise more easily
    #mol = doCanonicalAtomRenumbering(mol)
    getMolImageWithAtomNumbers(mol, imgFileName)
    
    if returnSmiles:
        #Choose atom idx
        atom_idx = int(input(f'Please select an atom index for the fragment exit vector: An image showing the index associated with each atom can be found at {imgFileName}.'))
        molWithDummy = addDummyAtomToMol(mol, atom_idx)
        removeAtomIndices(molWithDummy)
        smiles = Chem.MolToSmiles(molWithDummy)
        
        smiles = reformatDummyAtomInSMILES(smiles)
        
        if smilesFileName is not None:
            with open(smilesFileName, 'w') as f:
                f.write(smiles)
        
        print(f'You can use the following SMILES string as an input to STRIFE: {smiles}')

    else:
        print(f'An image showing the index associated with each atom can be found at {imgFileName}.')
        
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
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--fragment_sdf', '-f', type = str, required = True,
                        help = 'Location of fragment SDF.')
    parser.add_argument('--output_image_file', '-o', type = str, required = True,
                        help = 'Location to save numbered image of molecule')
    
    
    parser.add_argument('--return_smiles_string', '-r', action = "store_true", 
                        help = 'Return a SMILES string to be provided to STRIFE. If you do not request a SMILES string then you should provide the desired atom index as an argument to STRIFE.')    
    parser.add_argument('--save_smiles_string', '-s', default = None,
                        help = 'Location to save returned STRING (optional, SMILES string will be printed to the console)')
    
    
    arguments = parser.parse_args()
    
    molSDF = arguments.fragment_sdf
    imgFileName = arguments.output_image_file
    returnSmiles = arguments.return_smiles_string
    smilesFileName = arguments.save_smiles_string
    

    getSMILESStringForSTRIFE(molSDF, imgFileName, returnSmiles, smilesFileName)
    
    
    
    


    
    
    
    
    

