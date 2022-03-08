#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 09:42:58 2021

@author: Tom Hadfield

Class for generating molecules with STRIFE

"""


#Import Libaries

import warnings
warnings.filterwarnings("ignore",category=DeprecationWarning)


#########Standard Libraries##########
import json
import numpy as np
import pandas as pd
import tempfile
import os


#########RDKit Modules############
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') #Supress annoying RDKit output


from fineGrainedGenModel import DenseGGNNChemModel as DenseGGNNChemModelPharm
from GenModelDeLinker import DenseGGNNChemModel as DenseGGNNChemModelOrig
from coarseGrainedGenModel import DenseGGNNChemModel as DenseGGNNChemModelCount
from utils import elab_utils as eUtils


#TODO - review use of these
from utils import graph_utils, add_pharm_profile, frag_utils

import utils.single_frag_SMI2JSON_Pharm as s2JSON

#import summaryStats


class elaborate:
    
    #Collection of methods to support elaborations
    def __init__(self, useDefaultModels = True, models = None):
        
        #Specify the file names of the saved models to be used
        if useDefaultModels:
            self.savedModelDict = {'Orig':'models/GenModelDeLinker_saved.pickle',
                 'Count':'models/coarseGrainedGenModel_saved.pickle',
                 'Pharm':'models/fineGrainedGenModel_saved.pickle'}
        else:
            self.savedModelDict = models
            
    def addAlkane(self, fragment, alkLength):

        #method which allows us to generate elaborations of a pre-specified length using the DeLinker architecture
        #Under the original DeLinker architecture a seed molecule is used to derive attributes about the kinds of molecules
        #we want to elaborate (e.g. the elaboration length, desired pharmacophoric profile). In STRIFE, we specify those 
        #attributes ourselves without a seed molecule, but we use the existing DeLinker architecture to generate elaborations
        #of that length, by using a molecule comprising of the fragment and an alkane of the desired length as the seed molecule.
        
        
        #Create alkane to add to the fragment
        alk = 'C'*int(alkLength) + '[*:1]'
    
        fragAlk = f'{fragment}.{alk}'
        Chem.MolFromSmiles(fragAlk)
    
        #Identify dummy atoms and their neighbors
        m = Chem.MolFromSmiles(fragAlk)
        dummies = []
        dummyNeighbors = []
        for atom in m.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummies.append(atom.GetIdx())
                dummyNeighbors.append([n.GetIdx() for n in atom.GetNeighbors()][0])
    
    
    
        #Break the bonds attaching the dummy atoms and add a bond between their neighbors
        mw = Chem.RWMol(m)
        mw.RemoveBond(dummies[0], dummyNeighbors[0])
        mw.RemoveBond(dummies[1], dummyNeighbors[1])
    
        mw.AddBond(dummyNeighbors[0], dummyNeighbors[1], Chem.BondType.SINGLE)
    
        Chem.SanitizeMol(mw)
    
        #Now we just want to keep the full molecule and get rid of the dummies:
        smiles = Chem.MolToSmiles(mw).split('.')
        for s in smiles:
            if Chem.MolFromSmiles(s).GetNumHeavyAtoms() > 1:
                fullAlk = s
    
        return fullAlk
        
        
    def createSmi(self, fragment, alkLength):

        full = self.addAlkane(fragment, alkLength)
        chopped = 'C'*int(alkLength) + '[*:1]'
        smi = pd.DataFrame({'full':full, 'chopped':chopped, 'frag':fragment, 'dist':[0], 'ang':[0]})
    
        return smi
    
    def preprocessSingle(self, raw_data):
        processed_data = []
        dataset = 'zinc'
        
        for i, (smiles_mol, smiles_frag, QED_in, abs_dist) in enumerate([(mol['smiles_mol'], mol['smiles_frag'],
                                                                               mol['QED_in'], mol['abs_dist']) for mol in raw_data]):
           
            (mol_out, mol_in), nodes_to_keep, exit_points = eUtils.align_smiles_by_frags(smiles_mol, smiles_frag)
            if mol_out == []:
                continue
            nodes_in, edges_in = graph_utils.to_graph_mol(mol_in, dataset)
            nodes_out, edges_out = graph_utils.to_graph_mol(mol_out, dataset)
            if min(len(edges_in), len(edges_out)) <= 0:
                continue
            processed_data.append({
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
    
        return processed_data
    
    def replacePharmProfile(self, newProfile, jsonF, replacementType='Orig'):

        if replacementType not in ['Orig', 'Count', 'Pharm']:
            print('New pharmacophoric profile must be either Orig/Count/Pharm')
            return -1
    
        if replacementType == 'Orig' and len(newProfile) != 2:
            print('Profile must be [0, 0]')
            return -2
    
        if replacementType == 'Count' and len(newProfile) != 5:
            print('Profile must be of length 5')
            return -3
    
        if replacementType == 'Pharm' and len(newProfile) != 26:
            print('Profile must be of length 26')
            return -4
    
        #Replace profile
        jsonF[0]['abs_dist'] = newProfile

        return jsonF
    
    
    def generateMolsFromSpecificProfile(self, frag, elabLength, modelType, profile, numToGenerate = 250, smiName = 'outSmi'):

    #First want to create a model outside of the function so that we only need to update the model.valid_data attribute
    #Instead of loading the entire model repeatedly

        tempDir = tempfile.TemporaryDirectory()


        smi = self.createSmi(frag, elabLength)
        smi.to_csv(f'{tempDir.name}/{smiName}.smi', header = False, index = False, sep = ' ')
    
        smi.columns = ['full', 'chopped', 'frag', 'dist', 'ang']
        smi['pharm'] = [add_pharm_profile.createPharmacophoricProfile(row['frag'], row['full']) for idx, row in smi.iterrows()]
        smi.to_csv(f'{tempDir.name}/{smiName}Pharm.smi', header = None, index = None, sep = ' ')
    
        raw_data = s2JSON.train_valid_split(f'{tempDir.name}/{smiName}Pharm.smi')['valid']
    
        jsonFile = self.preprocessSingle(raw_data)
    
        #Replace pharmacophoric profile
        jsonFile = self.replacePharmProfile(profile, jsonFile, modelType)
    
        jsonName = f'{tempDir.name}/{smiName}Json{modelType}.json'
    
    
        with open(jsonName, 'w') as f:
            json.dump(jsonFile, f)
    
        data_input = jsonName
        
        config_input = '{"generation":true, "number_of_generation_per_valid": ' + str(numToGenerate) + ', "num_timesteps": 9, "hidden_size": 50, "encoding_size": 4, "batch_size": 1, "qed_trade_off_lambda": 0, "node_keep_trade_off_lambda": 0, "fwd_bkd_trade_off_lambda": 0, "num_epochs": 2, "epoch_to_generate":2, "compensate_num": 0, "num_different_starting": 1, "train_file":' + '"' + data_input + '"' + ', "valid_file":' + '"' + data_input + '"' + ', "output_dir": "' + tempDir.name + '", "random_seed": 4}'

        args = {'--dataset':'zinc', '--config':config_input,  '--restore': self.savedModelDict[modelType]}
    
    
    
        if modelType == 'Orig':
            model = DenseGGNNChemModelOrig(args)
        elif modelType == 'Count':
            model = DenseGGNNChemModelCount(args)
        elif modelType == 'Pharm':
            model = DenseGGNNChemModelPharm(args)
    

        #Generate Molecules
        model.train()
    
        #Get elaborations into a dataframe
        outDF = pd.read_csv(f"{tempDir.name}/{model.run_id}_generated_smiles_{model.params['dataset']}.smi", header = None, sep = ' ')
        outDF.columns = ['frag', 'full', 'gen']
        
        #Remove tempDir
        tempDir.cleanup()
        
        #print(outDF)

        return outDF
    
    
    def getPharmNumbers(self, mol):
    
        pharms = ['Acceptor', 'Donor', 'Aromatic']
        pharmCounts = {'Acceptor': 0, 'Donor': 0, 'Aromatic': 0}
        
    
        fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
        factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
    
        feats = factory.GetFeaturesForMol(mol)
        proc_feats = []
        for feat in feats:
            if feat.GetFamily() in pharms:
                proc_feats.append(feat.GetFamily())
                
        for p in proc_feats:        
            pharmCounts[p] += 1
        
        return pharmCounts
    
    
    def checkProfileCount(self, mol, frag, profile):
    #Take a molecule as input and check whether it has the same number of pharmacophores as specified in the profile

    #Profile should be of the form: [dist (=0), ang (=0), HBANum, HBDNum, AromaticNum]
    
        molPharms = self.getPharmNumbers(mol)
        fragPharms = self.getPharmNumbers(frag)
        
        elabProfile = [0, 0, molPharms['Acceptor'] - fragPharms['Acceptor'], molPharms['Donor'] - fragPharms['Donor'], molPharms['Aromatic'] - fragPharms['Aromatic']]
        
        if elabProfile == profile:
            return True
        else:
            return False
    

    def filterGeneratedMols(self, df, n_cores = 1):

        df = eUtils.parallelize_dataframe(df, eUtils.addElab, n_cores = n_cores)
        df = frag_utils.check_2d_filters_dataset_pandas(df)
    
        df = df.loc[df['passed'] == 1][['frag', 'full', 'gen', 'elab']]
    
        return df
    

    def makeElaborationsAndFilter(self, Hotspot, jobName = 'target', numElabsPerPoint = 250, n_cores = 1):

        elabLengths, profiles = Hotspot.determineProfiles()

        Hotspot.profileElabs = pd.DataFrame()
        for idx in range(len(elabLengths)):
            
            #get smiles of from file
            with open(Hotspot.frag, 'r') as f:
                fragSmiles = f.read()
                
            print(elabLengths[idx], profiles[idx], fragSmiles)

            g = self.generateMolsFromSpecificProfile(fragSmiles, elabLengths[idx], 'Count', profiles[idx], numToGenerate = numElabsPerPoint, smiName=jobName)
            #Filter:
            g = self.filterGeneratedMols(g, n_cores = n_cores)

            #filter those that don't have the correct profile out
            g['desiredProfile'] = [self.checkProfileCount(Chem.MolFromSmiles(row['gen']), Chem.MolFromSmiles(row['frag']), profiles[idx]) for i, row in g.iterrows()]
            g = g.loc[g['desiredProfile']]

            #Append g onto the list of molecules to be docked:
            Hotspot.profileElabs = Hotspot.profileElabs.append(g[['frag', 'full', 'gen']])

        return Hotspot
    
    
    def makeElaborationsNoFilter(self, Hotspot, modelType = 'Count', jobName = 'target', numElabsPerPoint = 250, filterQuality = False, n_cores = 1):
        
        #if filterQuality, then we filter the generated molecules using the set of 2D filters but we don't ensure they have the correct pharmacophoric profile
        
        elabLengths, profiles = Hotspot.determineProfiles()
        
        
        
        Hotspot.profileElabs = pd.DataFrame()
        for idx in range(len(elabLengths)):
            
            #get smiles of from file
            with open(Hotspot.frag, 'r') as f:
                fragSmiles = f.read()
            
            if modelType == 'Orig':
                profiles[idx] = [0, 0]
                
            print(elabLengths[idx], profiles[idx], fragSmiles)
            
        
        
            g = self.generateMolsFromSpecificProfile(fragSmiles, elabLengths[idx], modelType, profiles[idx], numToGenerate = numElabsPerPoint, smiName=jobName)
            print(g.shape[0])

            if filterQuality:
                g = self.filterGeneratedMols(g, n_cores = n_cores)
            print(g.shape[0])
            #Append g onto the list of molecules to be docked:
            Hotspot.profileElabs = Hotspot.profileElabs.append(g[['frag', 'full', 'gen']])

        return Hotspot
    
    
    
    def makePharmElabsQuasiActives(self, quasiActives, frag, filterMols = True, n_cores = 1):
        #make elaborations using the quasi actives extracted from a single hotspot or hotspot pair


        elabLengths = []
        pharmProfiles = []
        
        for s in quasiActives['smiles']:
            
            el = Chem.MolFromSmiles(s).GetNumHeavyAtoms() - Chem.MolFromSmiles(frag).GetNumHeavyAtoms()
            if el > 1:
                #Can't currently specify the desired elaboration length to just be one atom:
                elabLengths.append(el)
                p = add_pharm_profile.createPharmacophoricProfile(frag, s)
                p = add_pharm_profile.pharmProfileToList(p)
                pharmProfiles.append(p)
           

        pharmElabs = pd.DataFrame()

        for idx, p in enumerate(pharmProfiles):
            d = self.generateMolsFromSpecificProfile(frag, elabLengths[idx], 'Pharm', p, smiName = 'pharmElabs')
            
            if filterMols:
                #Apply 2D Filter:
                d = self.filterGeneratedMols(d, n_cores = n_cores)
            pharmElabs = pharmElabs.append(d)
            
        return pharmElabs
