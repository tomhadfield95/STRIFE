#Takes a smi file as input. Adds a sixth column on the end where we include the pharmacophoric profile


#Import data and modules
import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit.Chem import rdMolAlign
from rdkit.Chem import QED
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMMPA
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit import DataStructs
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolAlign

import matplotlib.pyplot as plt
from numpy import linalg
from scipy import stats
from multiprocessing import Pool
from random import sample


import os
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
import sys
import json
from random import sample

from rdkit import DataStructs
#from ScoringFunctionsTEMPORARY import *



#sys.argv[1] a smi file (without the extension)


#Create feature factory
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)



def align_smiles_by_frags(smiles_mol, smiles_frag):
    #Amended function which takes a single fragment as input
    try:
        smiles_frags = smiles_frag + '.[*:2]'
        mols_to_align = [Chem.MolFromSmiles(smiles_mol), Chem.MolFromSmiles(smiles_frags)]
        frags = [Chem.MolFromSmiles(smiles_frag)]

        # Include dummy in query
        du = Chem.MolFromSmiles('*')
        qp = Chem.AdjustQueryParameters()
        qp.makeDummiesQueries=True

        # Renumber based on frags (incl. dummy atoms)
        aligned_mols = []
        for i, mol in enumerate(mols_to_align):
            sub_idx = []
            for frag in frags:
                # Align to frags
                qfrag = Chem.AdjustQueryProperties(frag,qp)
                sub_idx += list(mol.GetSubstructMatch(qfrag))
            nodes_to_keep = [i for i in range(len(sub_idx))]
            if i == 0:
                mol_range = list(range(mol.GetNumHeavyAtoms()))
            else:
                mol_range = list(range(mol.GetNumHeavyAtoms()+2))
            idx_to_add = list(set(mol_range).difference(set(sub_idx)))
            sub_idx.extend(idx_to_add)
    #         print(sub_idx)
    #         print(mol.GetNumAtoms())
            aligned_mols.append(Chem.rdmolops.RenumberAtoms(mol, sub_idx))

        # Renumber dummy atoms to end
        dummy_idx = []
        for atom in aligned_mols[1].GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummy_idx.append(atom.GetIdx())
        #print(dummy_idx)
        for i, mol in enumerate(aligned_mols):
            sub_idx = list(range(aligned_mols[1].GetNumHeavyAtoms()+2))
            for idx in dummy_idx:
                sub_idx.remove(idx)
                sub_idx.append(idx)
            if i == 0:
                mol_range = list(range(mol.GetNumHeavyAtoms()))
            else:
                mol_range = list(range(mol.GetNumHeavyAtoms()+2))
            idx_to_add = list(set(mol_range).difference(set(sub_idx)))
            sub_idx.extend(idx_to_add)
            #print(sub_idx)
            #print(mol.GetNumAtoms())
            aligned_mols[i] = Chem.rdmolops.RenumberAtoms(mol, sub_idx[0:mol.GetNumAtoms()])

            # Get exit vectors
        exit_vectors = []
        for atom in aligned_mols[1].GetAtoms():
            if atom.GetAtomicNum() == 0:
                if atom.GetIdx() in nodes_to_keep:
                    nodes_to_keep.remove(atom.GetIdx())
                for nei in atom.GetNeighbors():
                    exit_vectors.append(nei.GetIdx())

        if len(exit_vectors) != 1:
            print("Incorrect number of exit vectors")

        return (aligned_mols[0], aligned_mols[1]), nodes_to_keep, exit_vectors

    except:
        print("Could not align")
        return ([],[]), [], []

def remove_dummy_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol2 = AllChem.ReplaceSubstructs(mol, Chem.MolFromSmiles('*'), Chem.MolFromSmiles('[H]'), True)[0]
    mol3 = Chem.RemoveHs(mol2)

    return Chem.MolToSmiles(mol3)


def HBondDonor(genSmiles, fragSmiles, factory = factory):
    #Input two smiles strings,
    #Align molecules and check whether a HBD exists in the part 
    #of the supermolecule which isn't in the fragment

    (aligned_mol, aligned_frag), to_keep, exit = align_smiles_by_frags(genSmiles, fragSmiles)

    feats = factory.GetFeaturesForMol(aligned_mol) #Gets all pharmacophoric features
    donors_idx = []

    for feat in feats:
        if feat.GetType() == 'SingleAtomDonor' and feat.GetAtomIds()[0] not in to_keep:
            donors_idx.append(feat.GetAtomIds()[0])

    path_lengths = []
    for donor in donors_idx:
        path = rdmolops.GetShortestPath(aligned_mol, exit[0], donor)
        path_lengths.append(len(path))

    #path_lengths contains the distance of each HBD to the exit point
    #If it has length 0 then we know there are no HBD in the elaborated structure
    return path_lengths

def HBondAcceptor(genSmiles, fragSmiles, factory = factory):
    #Input two smiles strings,
    #Align molecules and check whether a HBD exists in the part
    #of the supermolecule which isn't in the fragment

    (aligned_mol, aligned_frag), to_keep, exit = align_smiles_by_frags(genSmiles, fragSmiles)

    feats = factory.GetFeaturesForMol(aligned_mol) #Gets all pharmacophoric features
    acceptors_idx = []

    for feat in feats:
        if feat.GetType() == 'SingleAtomAcceptor' and feat.GetAtomIds()[0] not in to_keep:
            acceptors_idx.append(feat.GetAtomIds()[0])

    path_lengths = []
    for acceptor in acceptors_idx:
        path = rdmolops.GetShortestPath(aligned_mol, exit[0], acceptor)
        path_lengths.append(len(path))

    #path_lengths contains the distance of each HBD to the exit point
    #If it has length 0 then we know there are no HBD in the elaborated structure
    return path_lengths


def locAromatic(genSmiles, fragSmiles):
    #Function which checks whether any of the elaborated atoms belong in an aromatic system.
    #Returns a list containing 0 if there are no aromatic groups in the elaborated structure
    #Returns a list containing the minimal path length if there is an aromatic group
    (aligned_mol, aligned_frag), to_keep, exit = align_smiles_by_frags(genSmiles, fragSmiles)


    aromatic_idx = []

    for atom in aligned_mol.GetAtoms():
        if atom.GetIdx() not in to_keep and aligned_mol.GetAtomWithIdx(atom.GetIdx()).GetIsAromatic():
            #If the atom is in the elaborated structure and is aromatic
            aromatic_idx.append(atom.GetIdx())

    path_lengths = []
    for arom in aromatic_idx:
        path_lengths.append(len(rdmolops.GetShortestPath(aligned_mol, exit[0], arom)))

    if len(aromatic_idx) > 0:
        return [min(path_lengths)]
    else:
        return [0]


def createPharmacophoricProfile(frag, full):
    
    outDict = {'numHBA':0, 'numHBD':0, 'Aromatic':0, 'HBAPath':[0]*10, 'HBDPath':[0]*10, 'AroPath':[0]}
    
    try:

        #Get Number and location of HBA
        HBAPath = HBondAcceptor(full, frag)
        outDict['numHBA'] = len(HBAPath)
        
        if len(HBAPath) > 0:
            for idx, p in enumerate(HBAPath):
                outDict['HBAPath'][idx] = p
        

        
        #Get Number and location of HBD
        HBDPath = HBondDonor(full, frag)
        outDict['numHBD'] = len(HBDPath)
        
        if len(HBDPath) > 0:
            for idx, p in enumerate(HBDPath):
                outDict['HBDPath'][idx] = p
        
        
        #Ascertain presence and location of Aromatic group:
        arom = locAromatic(full, frag)
        
        if arom[0] > 0: #i.e. there is an aromatic group
            outDict['Aromatic'] = 1
            outDict['AroPath'] = arom
            
        outList = [outDict['numHBA'], outDict['numHBD'], outDict['Aromatic'], outDict['HBAPath'], outDict['HBDPath'], outDict['AroPath']]
        
        return outList

    except:

        return -1

def pharmProfileToList(profile):
    #Take as input a pharmacophoric profile of the form [NumHBA, NumHBD, NumAromatic, [HBAPath], [HBDPath], [Aromatic Path]]
    #And return a single list (i.e. take out the sublists)
    
    return [0,0] + profile[0:3] + profile[3] + profile[4] + profile[5] #[0,0] includes dist and ang placeholders


if __name__ == "__main__":

    smiData = pd.read_csv(sys.argv[1] + '.smi', header = None, sep = ' ')
    smiData.columns = ['full', 'chopped', 'frag', 'dist', 'ang']
    
    #Add Pharmacophoric features to the smi data
    smiData['pharm'] = [createPharmacophoricProfile(row['frag'], row['full']) for idx, row in smiData.iterrows()]
    
    smiData = smiData.loc[smiData['pharm'] != -1]

    #Write to file
    smiData.to_csv(sys.argv[1] + 'Pharm.smi', header = False, index = False, sep = ' ')




