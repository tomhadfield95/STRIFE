#!/usr/bin/env/python
"""
Usage:
    get_paired_zinc_test.py

Options:
    -h --help                Show this screen.


This script should take two sys.argv arguments
1) A SMI file 
2) A string which will be used as the root for the names of the output

Alongside angle and distance, we add our pharmacophoric profile
so the angle and distance vector will now bw:

[distance, angle, numHBA, numHBD, numAromatic, location of HBA, location of HBD, location of first aromatic atom]
The location of HBD is a list of 10 numbers (padded with zeros) but we we don't pass it here as a list (i.e. we're still working with a one dimensional object)
"""

import sys, os

#COMMENTING THIS LINE OUT FOR USE IN noActives pipeline
#If looking to use this script for standard fragment elaboration then need to uncomment this line
#sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import QED
from rdkit.Chem import rdFMCS
import glob
import json
import numpy as np
from utils import bond_dict, dataset_info, need_kekulize, to_graph, to_graph_mol, align_smiles_by_MCS, graph_to_adj_mat
import utils
import pickle
import random
from docopt import docopt
from align_molecules import align_smiles_by_MCS_repeat
import pandas as pd

dataset = 'zinc'

def read_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    num_lines = len(lines)
    data = []
    for i, line in enumerate(lines):
        smiles_mol, smiles_rest, smiles_frag, abs_dist, angle, Pharm = line.strip().split(' ')
        #smiles_mol, smiles_link, smiles_frag1, smiles_frag2, abs_dist, angle, _, _ = line.strip().split(' ')
        mu_mol = QED.qed(Chem.MolFromSmiles(smiles_mol))

        numHBA = Pharm[0]
        numHBD = Pharm[1]
        numAromatic = Pharm[2]
        pathHBA = Pharm[3]
        pathHBD = Pharm[4]
        pathAromatic = Pharm[5]


        data.append({'smiles_mol': smiles_mol, 'smiles_rest' : smiles_rest, 
                     'smiles_frag': smiles_frag,
                     'QED_in': mu_mol, 'abs_dist': [abs_dist,angle, numHBA, numHBD, numAromatic] + pathHBA + pathHBD + pathAromatic})
        if i % 2000 ==0:
            print('finished reading: %d / %d' % (i, num_lines), end='\r')
    return data

def read_file_pandas(file_path):
    smi_data = pd.read_csv(file_path, header = None, sep = " ")
    smi_data.columns = ['full', 'chopped', 'frag', 'dist', 'ang', 'Pharm']
    data = []

    for idx, row in smi_data.iterrows():
        
        mu_mol = QED.qed(Chem.MolFromSmiles(row['full']))
        Pharms = json.loads(row['Pharm'])
        data.append({'smiles_mol': row['full'], 'smiles_rest': row['chopped'], 'smiles_frag': row['frag'], 'QED_in':mu_mol, 'abs_dist': [row['dist'], row['ang'], Pharms[0], Pharms[1], Pharms[2]] + Pharms[3] + Pharms[4] + Pharms[5]})

    return data

'''
def read_file_test(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    num_lines = len(lines)
    data = []
    for i, line in enumerate(lines):
        smiles_in = line.strip().split(' ')[0]
        mu_in = QED.qed(Chem.MolFromSmiles(smiles_in))
        data.append({'smiles_in': smiles_in, "smiles_out": smiles_in,
                     'QED_in': mu_in, 'QED_out': mu_in})
        if i % 2000 ==0:
            print('finished reading: %d / %d' % (i, num_lines), end='\r')
    return data
'''

def train_valid_split(smiles_path):
    print('reading data...')
    raw_data = {'train': [], 'valid': []} # save the train, valid dataset.

    all_data = read_file_pandas(smiles_path)

    size_data = len(all_data)
    #valid_idx = int(size_data*0.95) # 90/10 train/valid split. Sequential split.
    valid_idx = len(all_data) - 1000 # Fixed 1000 examples for validation
    #valid_idx = 10000 # arbitrary large number

    file_count=0
    # Below inefficient given split type, but leaving for if changed
    for i, data_item in enumerate(all_data):
        if i < valid_idx:
            raw_data['train'].append(data_item)
        else:
            raw_data['valid'].append(data_item)
    return raw_data

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
            #print(sub_idx)
            #print(mol.GetNumHeavyAtoms())
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
            aligned_mols[i] = Chem.rdmolops.RenumberAtoms(mol, sub_idx)

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


'''
def align_smiles_by_frags(smiles_mol, smiles_frag1, smiles_frag2):
    try:
        smiles_frags = smiles_frag1 + '.' + smiles_frag2
        mols_to_align = [Chem.MolFromSmiles(smiles_mol), Chem.MolFromSmiles(smiles_frags)]
        frags = [Chem.MolFromSmiles(smiles_frag1), Chem.MolFromSmiles(smiles_frag2)]

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
            #print(sub_idx)
            #print(mol.GetNumHeavyAtoms())
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
            aligned_mols[i] = Chem.rdmolops.RenumberAtoms(mol, sub_idx)

        # Get exit vectors
        exit_vectors = []
        for atom in aligned_mols[1].GetAtoms():
            if atom.GetAtomicNum() == 0:
                if atom.GetIdx() in nodes_to_keep:
                    nodes_to_keep.remove(atom.GetIdx())
                for nei in atom.GetNeighbors():
                    exit_vectors.append(nei.GetIdx())

        if len(exit_vectors) != 2:
            print("Incorrect number of exit vectors")

        return (aligned_mols[0], aligned_mols[1]), nodes_to_keep, exit_vectors

    except:
        print("Could not align")
        return ([],[]), [], []
'''

def preprocess(raw_data, dataset, name):
    print('\nparsing smiles as graphs...')
    processed_data = {'train': [], 'valid': []}
    
    file_count = 0
    for section in ['train', 'valid']:
        all_smiles = [] # record all smiles in training dataset
        for i, (smiles_mol, smiles_frag, QED_in, abs_dist) in enumerate([(mol['smiles_mol'], mol['smiles_frag'], 
                                                                               mol['QED_in'], mol['abs_dist']) for mol in raw_data[section]]):
            #print(file_count, end='\r')
            #print("File %d" % i)
            (mol_out, mol_in), nodes_to_keep, exit_points = align_smiles_by_frags(smiles_mol, smiles_frag)
            if mol_out == []:
                continue
            nodes_in, edges_in = to_graph_mol(mol_in, dataset)
            nodes_out, edges_out = to_graph_mol(mol_out, dataset)
            if min(len(edges_in), len(edges_out)) <= 0:
                continue
            processed_data[section].append({
                'targets_in': [[(QED_in)]],
                'targets_out': [[(QED_in)]],
                'graph_in': edges_in,
                'graph_out': edges_out, # in to ensure same number of v in in and out - HACK
                'node_features_in': nodes_in,
                'node_features_out': nodes_out, # in to ensure same number of v in in and out - HACK
                'smiles_out': smiles_mol,
                'smiles_in': smiles_frag + '.[*:2]',
                'v_to_keep': nodes_to_keep,
                'exit_points': exit_points,
                'abs_dist': abs_dist
            })
            all_smiles.append(smiles_mol)
            all_smiles.append(smiles_frag + '.[*:2]')
            if file_count % 500 == 0:
                print('finished processing: %d' % file_count, end='\r')
            file_count += 1
        print('\n%s: 100 %%      ' % (section))
        # save the dataset
        with open('molecules_%s_%s_%s.json' % (name, section, dataset), 'w') as f:
            print("Section length %d" % len(processed_data[section]))
            json.dump(processed_data[section], f)
        # save all molecules in the training dataset
        #if section == 'train':
        #    utils.dump('smiles_%s.pkl' % dataset, all_smiles)         
          
'''
def preprocess_test(smiles_path, dataset):
    print('loading data')
    raw_data = read_file_test(smiles_path)
    print('parsing smiles as graphs...')
    processed_data = []

    file_count = 0
    all_smiles = [] # record all smiles in training dataset
    for i,(smiles_in, smiles_out, QED_in, QED_out) in enumerate([(mol['smiles_in'], mol['smiles_out'],
                                        mol['QED_in'], mol['QED_out']) for mol in raw_data]):
        (mol_in, mol_out), _, nodes_to_keep = align_smiles_by_MCS(smiles_in, smiles_out) # Fine as mols identical
        nodes_in, edges_in = to_graph_mol(mol_in, dataset)
        nodes_out, edges_out = to_graph_mol(mol_out, dataset)
        if min(len(edges_in), len(edges_out)) <= 0:
            continue
        processed_data.append({
            'targets_in': [[(QED_in)]],
            'targets_out': [[(QED_out)]],
            'graph_in': edges_in,
            'graph_out': edges_out,
            'node_features_in': nodes_in,
            'node_features_out': nodes_out,
            'smiles_in': smiles_in,
            'smiles_out': smiles_out,
            'v_to_keep': nodes_to_keep
        })
        all_smiles.append(smiles_in)
        all_smiles.append(smiles_out)
        if file_count % 100 == 0:
            print('finished processing: %d' % file_count, end='\r')
        file_count += 1
    print('\n%test: 100 %%      ' )
    # save the dataset
    with open('molecules_%s_%s_%s.json' % ("nodepred", "test", dataset), 'w') as f:
        json.dump(processed_data, f)
'''

if __name__ == "__main__":
    data_path = sys.argv[1]
    name = sys.argv[2]
    #data_path = "fragments_casf_dist_ang.smi"
    #name = "fragments_casf_dist_ang_exact.smi"
    raw_data = train_valid_split(data_path)
    preprocess(raw_data, dataset, name)
