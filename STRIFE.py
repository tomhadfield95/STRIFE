#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 17:12:05 2021

@author: Tom Hadfield

Main class for STRIFE algorithm
"""



#########Standard Libraries##########
import json
import time
import argparse
import os

#import warnings
#warnings.filterwarnings("ignore",category=DeprecationWarning)
#warnings.filterwarnings("ignore",category=DeprecationWarning, module = 'tensorflow')

from elaborations import elaborate
import numpy as np
import pandas as pd
import multiprocessing as mp #For parallelising
import time
import matplotlib.pyplot as plt
import sys 
from random import sample
from functools import partial
import glob
import openbabel
from datetime import datetime
from data_prep.specifyExitVector import addDummyAtomToMol

#########RDKit Modules############
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import QED
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMMPA
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole 
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') #Supress annoying RDKit output

import pickle



from docking import docking

try:
    #Try to import the preprocessing class which uses the CSD package
    from preprocessing import preprocessing
    from hotspots.hs_io import HotspotReader
except:
    #If we can't do that import the class without those methods
    from preprocessing_no_hotspots import preprocessing



class STRIFE:
    
    
    def __init__(self, args):
        #run 'python STRIFE.py -h' for definitions of the arguments
        
        print('Running STRIFE Algorithm....')
       
        
        print('Doing argument checking...')
        
        assert bool(args.protein) == 1, 'Please specify the path to a PDB file'
        assert bool(args.fragment_sdf) == 1, 'Please specify the location of the fragment SDF. This can also be an SDF of a larger ligand of which the fragment is a substructure'
        assert bool(args.fragment_smiles) + bool(args.exit_vector_idx is not None) == 1, 'Please specify exactly one of: The location of a text file which contains the SMILES string of the fragment or a SMILES string (as the argument fragment_smiles) or the atomic index of the desired exit vector (as the argument exit_vector).'
        
        #Convert the provided paths to the absolute path 
        args.protein = os.path.abspath(os.path.expanduser(args.protein))
        args.fragment_sdf = os.path.abspath(os.path.expanduser(args.fragment_sdf))

        #Check that the arguments exist
        
        if args.protein is None:
            raise ValueError('You must specify a pdb file as the protein')
        else:
            assert os.path.exists(args.protein), f'Specified protein file, {args.protein}, does not exist'
        
        if args.fragment_sdf is None:
            raise ValueError('You must specify an SDF file, either containing the molecule to be used as a fragment, or a superstructure of it')
        else:
            assert os.path.exists(args.fragment_sdf), f'Specified fragment SDF file, {args.fragment_sdf}, does not exist'
        
        #If the output directory doesn't exist, create it
        if not os.path.exists(args.output_directory):
            os.makedirs(args.output_directory)

        args.output_directory = os.path.abspath(os.path.expanduser(args.output_directory))




        if args.fragment_smiles is not None:
            
            #Check whether args.fragment_smiles is a file
            if os.path.exists(os.path.expanduser(args.fragment_smiles)):
                args.fragment_smiles = os.path.abspath(os.path.expanduser(args.fragment_smiles))
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
                with open(args.fragment_smiles, 'r') as f:
                    fragSmiles = f.read()
            
            else:
                fragSmiles = args.fragment_smiles
            
        
        else:
            #using exit_vector_idx instead
            mol = Chem.SDMolSupplier(args.fragment_sdf)[0]
            molWithDummyAtom = addDummyAtomToMol(mol, args.exit_vector_idx)
            fragSmiles = Chem.MolToSmiles(molWithDummyAtom)
        
        
        
        if args.write_elaborations_dataset:
            self.writeFinalElabs = True
        else:
            self.writeFinalElabs = False
        
        
     
        assert args.model_type in [0, 1, 2], 'Please provide a valid setting for generating molecules: 0, 1 or 2. See the "Running STRIFE" section of the github readme for more details'
        
        
        #Check that we're inputting a valid pharmacophoric representation
        assert bool(args.hotspots_output) + bool(args.calculate_hotspots) + bool(args.load_specified_pharms) == 1, 'Please specify exactly one way for STRIFE to incorporate structural information. Either provide an already calculated FHM, request that STRIFE calculates an FHM, or provide your own pharmacophoric points'
        

        print('Argument checking complete.')
        
        print('Processing pharmacophoric information')
        
        
        self.storeLoc = args.output_directory
        
        if args.num_cpu_cores > 0:
            self.num_cpu_cores = args.num_cpu_cores
        elif args.num_cpu_cores == -1:
            self.num_cpu_cores = mp.cpu_count()
        else:
            raise ValueError("Please supply a valid number of cores to use, or specify num_cpu_cores as -1 to use all available cores")


        #Create subclasses
        self.elaborate = elaborate()
        self.docking = docking()
        self.preprocessing = preprocessing()
        
        
        
        if args.hotspots_output is not None:
            self.hotspotsLoc = args.hotspots_output
        
        if bool(args.calculate_hotspots):
            print('Calculating Fragment Hotspot Map for input protein...\n')
            print('This may take a few minutes...\n')
                
            self.preprocessing.calculateFHM(args.protein, args.calculated_hotspots)
            self.hotspotsLoc = f'{args.calculated_hotspots}/out.zip'
            print(f'FHM calculation complete. Output saved at {self.hotspotsLoc}')
        
        
        if args.load_specified_pharms:
            
            #Load the specified pharmacophoric points from args.output_directory
            if len(glob.glob(f'{self.storeLoc}/acceptorHotspot.sdf')) + len(glob.glob(f'{self.storeLoc}/donorHotspot.sdf')) == 0:
                print('If manually specifying pharmacophoric points, please provide at least one point\n')
                print(f'Donors should be saved in the file {self.storeLoc}/donorHotspot.sdf\n')
                print(f'Acceptors should be saved in the file {self.storeLoc}/acceptorHotspot.sdf\n')
                raise ValueError('If manually specifying pharmacophoric points, please provide at least one point')
                
            elif len(glob.glob(f'{self.storeLoc}/acceptorHotspot.sdf')) == 0:
                hotspotsDict = {}
                hotspotsDict['Acceptor'] = Chem.RWMol()
                hotspotsDict['Donor'] = Chem.MolFromMolFile(f'{self.storeLoc}/donorHotspot.sdf')
            elif len(glob.glob(f'{self.storeLoc}/donorHotspot.sdf')) == 0:
                hotspotsDict = {}
                hotspotsDict['Acceptor'] = Chem.MolFromMolFile(f'{self.storeLoc}/acceptorHotspot.sdf')
                hotspotsDict['Donor'] = Chem.RWMol()
            else:
                hotspotsDict = {}
                hotspotsDict['Acceptor'] = Chem.MolFromMolFile(f'{self.storeLoc}/acceptorHotspot.sdf')
                hotspotsDict['Donor'] = Chem.MolFromMolFile(f'{self.storeLoc}/donorHotspot.sdf')
        
        
        
        print('Preprocessing fragment')
        
        #Store fragment SMILES in the output directory:
        with open(f'{self.storeLoc}/frag_smiles.smi', 'w') as f:
            f.write(fragSmiles)
        
            
        fragMol3D = Chem.SDMolSupplier(args.fragment_sdf)[0]
            
        fc, evI, evp, fc2 = self.preprocessing.preprocessFragment(fragSmiles, fragMol3D)
        
        #Save fragment SDF
        Chem.MolToMolFile(fc, f'{self.storeLoc}/frag.sdf')
        
        #Save constraint SDF (will need to be converted to Mol2 using obabel)
        Chem.MolToMolFile(fc2, f'{self.storeLoc}/constraint.sdf')
        
        #Save fragment exit position
        np.savetxt(f'{self.storeLoc}/evp.txt', evp)
        
        
        #Convert the constraint.sdf file to constraint.mol2 (for constrained docking in GOLD)
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("sdf", "mol2")
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, f'{self.storeLoc}/constraint.sdf')
        obConversion.WriteFile(mol, f'{self.storeLoc}/constraint.mol2')
        
        
        
        if args.load_specified_pharms == False:
            #i.e. we haven't already specified the pharmacophoric points
            
            #Read in hotspots output
            reader = HotspotReader(self.hotspotsLoc)
            reader.__enter__()
            results = reader.read()
             
            hotspotsDict = self.preprocessing.processHotspots(results, evp, fc)
            Chem.MolToMolFile(hotspotsDict['Acceptor'], f'{self.storeLoc}/acceptorHotspot.sdf')
            Chem.MolToMolFile(hotspotsDict['Donor'], f'{self.storeLoc}/donorHotspot.sdf')




        
        self.exitVectorPos = evp
        self.frag = f'{self.storeLoc}/frag_smiles.smi'
        self.fragCore = fc
        self.constraintFile = f'{self.storeLoc}/constraint.mol2'
        self.cavityLigandFile = args.fragment_sdf #Used for docking in GOLD to define the binding pocket
        self.protein = args.protein
        
        
        
        
        #Set up the hotspotsDF
        self.HPositions = [] #list for hotspot positions
        self.HType = [] #list for hotspot type (donor/acceptor)
        self.Distances = [] #Distance from exit vector
        self.Angles = [] #Angle from exit vector
        self.origAtomIdx = []

        for pharm in ['Acceptor', 'Donor']:

            if hotspotsDict[pharm].GetNumHeavyAtoms() > 0: #Changed self.hotspotsDict to hotspotsDict here
                for atom in hotspotsDict[pharm].GetAtoms():
                    pos = np.array(hotspotsDict[pharm].GetConformer().GetAtomPosition(atom.GetIdx()))

                    self.HPositions.append(pos)
                    self.Distances.append(self.preprocessing.vectorDistance(pos, self.exitVectorPos))
                    self.Angles.append(self.preprocessing.vectorAngle(pos, self.exitVectorPos))
                    self.HType.append(pharm)
                    self.origAtomIdx.append(atom.GetIdx()) #Atom index so we can recover it if necessary

        self.HotspotsDF = pd.DataFrame({'distFromExit':self.Distances, 'angFromExit':self.Angles, 'position':self.HPositions, 'type':self.HType}).sort_values('distFromExit', ascending = True).reset_index(drop = True)


        self.hSingles = self.preprocessing.prepareProfiles(self.HotspotsDF)
        self.hMulti = self.preprocessing.prepareProfiles(self.HotspotsDF, single = False) #For satisfying multiple pharmacophoric point simultaneously

        #Import pathLength classifier
        with open(args.path_length_model, 'rb') as f:
            self.clf = pickle.load(f)
            
    
    def run(self, args):
        
        #Function which runs a version of the STRIFE algorithm depending on the specified model_type
        if args.model_type == 0:
            self.runWholePipeline(totalNumElabs = args.number_elaborations, numElabsPerPoint = args.number_elaborations_exploration)
        elif args.model_type == 1:
            self.runCustomPharms(numElabsRefinement = args.number_elaborations, numElabsExploration = args.number_elaborations_exploration)
        elif args.model_type == 2:
            self.elaborationsWithoutRefinement(counts = True, totalNumElabs = args.number_elaborations, numElabsPerPoint = args.number_elaborations_exploration)
        else:
            print('Please specify a valid model setting')
    
            
    def runWholePipeline(self, totalNumElabs = 250, numElabsPerPoint = 250):
        
        if self.HotspotsDF.shape[0] > 0:
            self.exploration(numElabsPerPoint = numElabsPerPoint)
            self.identifyQuasiActives()
            
            numQuasiActives = sum([self.singleQuasiActives[k].shape[0] for k in self.singleQuasiActives.keys()])
            if numQuasiActives > 0:
                self.refinement(totalNumElabs)
                self.status = 'No Issues'
            else:
                self.status = 'No quasi-actives identified'
        else: 
            self.status = 'No suitable hotspots identified'


    def manuallySpecifyPharmPoints(self):
        #TODO - Implement
        pass
    
    def exploration(self, numElabsPerPoint = 250, n_cores = None):
        
        #Generate elaborations using the count model
        #Dock them using the constrained docking functionality in GOLD
        #Measure the distance between the pharmacophoric point and a matching 
        #pharmacophore in the elaboration.
        
        if n_cores is None:
            n_cores = self.num_cpu_cores

        self.singleElabs = {}
        self.singleDocks = {}
        self.singleDistances = {}
        
        for k in self.hSingles.keys():
    
            #iterate over the pharmacophoric points
            
            #set up HotspotSingle class
            self.singleElabs[k] = HotspotSingle(self.hSingles[k], self.frag, self.clf, self.constraintFile, self.cavityLigandFile, self.protein)
            
            #Make elaborations using the counts model and filter to retain those with the desired pharmacophoric profile
            self.singleElabs[k] = self.elaborate.makeElaborationsAndFilter(self.singleElabs[k], numElabsPerPoint=numElabsPerPoint, n_cores = n_cores)
            
            #Prepare the filtered elaborations to be docked in GOLD
            self.singleElabs[k] = self.docking.prepareForDocking(self.singleElabs[k], self.fragCore, f'{self.storeLoc}/countsElabs{k}.sdf')
            
            #do docking
            self.singleDocks[k] = self.docking.dockLigandsMP(self.singleElabs[k].dockingFname, self.constraintFile, self.cavityLigandFile, self.protein, n_processes = n_cores) #Dock in parallel
            
            #Compute distance to pharmacophoric point
            self.singleDistances[k] = self.docking.assessAllDocks(self.singleDocks[k], self.hSingles[k], True)


    def identifyQuasiActives(self):
        
        self.singleQuasiActives = {}
        
        for k in self.hSingles.keys():
            #Select up to 5 molecules where the pharmaocphore is within 2A of the pharmacophoric point
            self.singleQuasiActives[k] = self.singleDistances[k].loc[self.singleDistances[k]['distance'] < 2].drop_duplicates('smiles').head(5)

        
    
    def refinement(self, totalNumElabs = 250, n_cores = None):
        #Use the quasi-actives to generate elaborations using the pharm model
        
        if n_cores is None:
            n_cores = self.num_cpu_cores

        if totalNumElabs is None:
            #i.e. we want all the elaborations derived from the quasi-actives
            self.singlePharmElabs = {}

            for k in self.hSingles.keys():
                self.singlePharmElabs[k] = pharmElabs() #create pharmElabs class to store outputs
                self.singlePharmElabs[k].profileElabs = self.elaborate.makePharmElabsQuasiActives(self.singleQuasiActives[k], self.frag) #make elaborations using the quasi actives profiles
                self.singlePharmElabs[k] = self.docking.prepareForDocking(self.singlePharmElabs[k], self.fragCore, f'{self.storeLoc}/pharmsElabs{k}.sdf')
        
                #Dock
                self.singlePharmElabs[k].docks, self.singlePharmElabs[k].fitnessScores = self.docking.dockLigandsMP(self.singlePharmElabs[k].dockingFname, self.constraintFile, self.cavityLigandFile, self.protein, returnFitnessScore = True, n_processes = n_cores)
    
                #Calculate Ligand Efficiency
                self.singlePharmElabs[k].ligEff = self.docking.ligandEfficiency(self.singlePharmElabs[k].docks, self.singlePharmElabs[k].fitnessScores)
                
        
        
        else:
            #i.e. we've specified a certain number of elaborations to make 
            
            #Iterate over all of the single hotspots, make elabs from their quasi actives, append them together and randomly sample totalNumElabs of them
            self.singlePharmElabs = {}
            self.pharmElabsTest = pd.DataFrame()
    
            with open(self.frag, 'r') as f:
                fragSmiles = f.read()
    
            for k in self.hSingles.keys():
                self.singlePharmElabs[k] = pharmElabs()
                self.singlePharmElabs[k].profileElabs = self.elaborate.makePharmElabsQuasiActives(self.singleQuasiActives[k], fragSmiles, filterMols = False) #make elaborations using the quasi actives profiles
                self.pharmElabsTest = self.pharmElabsTest.append(self.singlePharmElabs[k].profileElabs)
    
            #Now we have all of the molecules sampled from the quasi active profiles
            self.pharmElabsTestSample = self.pharmElabsTest.sample(n = totalNumElabs, random_state = 10) #Sample the number of elaborations we want
            self.pharmElabsTestFilter = self.elaborate.filterGeneratedMols(self.pharmElabsTestSample, n_cores = n_cores)
    
            #Now prepare for docking 
            self.pharmElabsTestFName = f'{self.storeLoc}/pharmsElabsTestPreDocking.sdf'
            self.pharmElabsCountsIdx, self.pharmElabsCountsSDFIdx = self.docking.getSDFs(self.pharmElabsTestFilter, self.fragCore, self.pharmElabsTestFName) 

            #Do docking
            self.pharmElabsTestDocks, self.pharmElabsTestFS = self.docking.dockLigandsMP(self.pharmElabsTestFName, self.constraintFile, self.cavityLigandFile, self.protein, returnFitnessScore = True, n_processes = n_cores)

            #Compute ligand Efficiency
            self.pharmElabsTestLigEff = self.docking.ligandEfficiency(self.pharmElabsTestDocks, self.pharmElabsTestFS)

            #Write the docks to file with the ligand efficiency as an attribute
            self.pharmElabsDockedFName = f'{self.storeLoc}/pharmsElabsTestDocked.sdf'
            w = Chem.SDWriter(self.pharmElabsDockedFName)
            
            for idx, m in enumerate(self.pharmElabsTestDocks):
                m.SetProp('STRIFE_LigEff_Score', str(self.pharmElabsTestFS[idx]/m.GetNumHeavyAtoms()))
                w.write(m)

            w.close()

            #Standardise the name of the final df
            self.rankedElaborationsFinal = self.pharmElabsTestLigEff

            if self.writeFinalElabs:
                self.rankedElaborationsFinal.to_csv(f'{self.storeLoc}/rankedElaborationsFinal.csv')

    
    def elaborationsWithoutRefinement(self, counts = True, totalNumElabs = 250, numElabsPerPoint = 250, n_cores = None):
        
        if n_cores is None:
            n_cores = self.num_cpu_cores

        if totalNumElabs == None:
            
            #Just
            
            self.singleElabs_noRefine = {}
            self.singleDocks_noRefine = {}
            
            for k in self.hSingles.keys():
                #iterate over the pharmacophoric points

                #set up HotspotSingle class
                self.singleElabs_noRefine[k] = HotspotSingle(self.hSingles[k], self.frag, self.clf, self.constraintFile, self.cavityLigandFile, self.protein)
                
                if counts:
                    #Make elaborations using the counts model
                    #Filter for elaboration quality before docking, but not on pharmacophoric profile
                    self.singleElabs_noRefine[k] = self.elaborate.makeElaborationsNoFilter(self.singleElabs_noRefine[k], numElabsPerPoint=numElabsPerPoint, filterQuality = True, n_cores = n_cores)
                else:
                    #Make elaborations using the Orig model 
                    self.singleElabs_noRefine[k] = self.elaborate.makeElaborationsNoFilter(self.singleElabs_noRefine[k], modelType = 'Orig', numElabsPerPoint=numElabsPerPoint, filterQuality = True, n_cores = n_cores)
                
                
                #Prepare the filtered elaborations to be docked in GOLD
                self.singleElabs_noRefine[k] = self.docking.prepareForDocking(self.singleElabs_noRefine[k], self.fragCore, f'{self.storeLoc}/countsElabsNoRefine{k}.sdf')
                
                
                #Now we want to dock and compute ligand efficiency
                self.singleElabs_noRefine[k].docks, self.singleElabs_noRefine[k].fitnessScores = self.docking.dockLigandsMP(self.singlePharmElabs[k].dockingFname, self.constraintFile, self.cavityLigandFile, self.protein, returnFitnessScore = True, n_processes = n_cores)
                    
                #Calculate Ligand Efficiency
                self.singleElabs_noRefine[k].ligEff = self.docking.ligandEfficiency(self.singleElabs_noRefine[k].docks, self.singleElabs_noRefine[k].fitnessScores) 
                
                
                #Write the docks to file with the ligand efficiency as an attribute
                w = Chem.SDWriter(f'{self.storeLoc}/countElabsNoRefine{k}_Docked.sdf')
                
                for idx, m in enumerate(self.singleElabs_noRefine[k].docks):
                    m.SetProp('STRIFE_LigEff_Score', str(self.singleElabs_noRefine[k].fitnessScores/m.GetNumHeavyAtoms()))
                    w.write(m)
    
                w.close()

            
            #Standardise the name of the final df
            self.rankedElaborationsFinal = pd.DataFrame()

            for k in self.hSingles.keys():
                self.rankedElaborationsFinal = self.rankedElaborationsFinal.append(self.singleElabs_noRefine[k].ligEff)
            
            if self.writeFinalElabs:
                self.rankedElaborationsFinal.to_csv(f'{self.storeLoc}/rankedElaborationsFinal.csv')


        else:
            
            self.singleElabs_noRefine = {}
            self.singleDocks_noRefine = {}

            for k in self.hSingles.keys():
                #iterate over the pharmacophoric points

                #set up HotspotSingle class
                self.singleElabs_noRefine[k] = HotspotSingle(self.hSingles[k], self.frag, self.clf, self.constraintFile, self.cavityLigandFile, self.protein)

                if counts:
                    #Make elaborations using the counts model
                    #Here we don't filter on quality until the final set of elaborations has been sampled
                    self.singleElabs_noRefine[k] = self.elaborate.makeElaborationsNoFilter(self.singleElabs_noRefine[k], numElabsPerPoint=numElabsPerPoint, n_cores = n_cores)
                else:
                    #Make elaborations using the Orig model 
                    self.singleElabs_noRefine[k] = self.elaborate.makeElaborationsNoFilter(self.singleElabs_noRefine[k], modelType = 'Orig', numElabsPerPoint=numElabsPerPoint, n_cores = n_cores)
        
    
            #Now we want to take all of the elaborations we've made, sample some of them and then filter
            self.elabsTestNoRefine = pd.DataFrame()
            for k in self.hSingles.keys():
                self.elabsTestNoRefine = self.elabsTestNoRefine.append(self.singleElabs_noRefine[k].profileElabs)
    
            #Now we have all of the molecules sampled from the quasi active profiles
            self.elabsTestNoRefineSample = self.elabsTestNoRefine.sample(n = totalNumElabs, random_state = 10) #Sample the number of elaborations we want
            self.elabsTestNoRefineFilter = self.elaborate.filterGeneratedMols(self.elabsTestNoRefineSample, n_cores = n_cores)
        
        
            #Now we want to dock these 
            self.elabsTestNoRefineFName = f'{self.storeLoc}/elabsTestNoRefine.sdf'
            self.elabsTestNoRefineCountsIdx, self.elabsTestNoRefineCountsSDFIdx = self.docking.getSDFs(self.elabsTestNoRefineFilter, self.fragCore, self.elabsTestNoRefineFName) 
            
            
            
            #Do docking
            self.elabsTestNoRefineDocks, self.elabsTestNoRefineFS = self.docking.dockLigandsMP(self.elabsTestNoRefineFName, self.constraintFile, self.cavityLigandFile, self.protein, returnFitnessScore = True, n_processes = n_cores)

            #Compute ligand Efficiency
            self.elabsTestNoRefineLigEff = self.docking.ligandEfficiency(self.elabsTestNoRefineDocks, self.elabsTestNoRefineFS)
    
    
            #Write the docks to file with the ligand efficiency as an attribute
            w = Chem.SDWriter(f'{self.storeLoc}/elabsTestNoRefine_Docked.sdf')
            
            for idx, m in enumerate(self.elabsTestNoRefineDocks):
                m.SetProp('STRIFE_LigEff_Score', str(self.elabsTestNoRefineFS/m.GetNumHeavyAtoms()))
                w.write(m)
    
                w.close()
    

            #Standardise the name of the final df
            self.rankedElaborationsFinal = self.elabsTestNoRefineLigEff

            if self.writeFinalElabs:
                self.rankedElaborationsFinal.to_csv(f'{self.storeLoc}/rankedElaborationsFinal.csv')
    
        
    def runCustomPharms(self, numElabsRefinement = 250, numElabsExploration = 250, n_cores = None):
        
        #Version of the STRIFE algorithm which supports running multiple pharmacophoric points simultaneously
        #It's recommended that you only use this with manually specified pharmacophoric points (or at least check those derived by the hotspots algorithm for suitability)
        
        if n_cores is None:
            n_cores = self.num_cpu_cores

        ########EXPLORATION PHASE##########
        self.multiElabs = HotspotMulti(self.hMulti, self.frag, self.clf, self.constraintFile, self.cavityLigandFile, self.protein)
        
        #Make elaborations using the counts model and filter to retain those with the desired pharmacophoric profile
        self.multiElabs = self.elaborate.makeElaborationsAndFilter(self.multiElabs, numElabsPerPoint=numElabsExploration, n_cores = n_cores)
        
        #Prepare the filtered elaborations to be docked in GOLD
        self.multiElabs = self.docking.prepareForDocking(self.multiElabs, self.fragCore, f'{self.storeLoc}/countsElabsMulti.sdf')
        
        #do docking
        self.multiDocks = self.docking.dockLigandsMP(self.multiElabs.dockingFname, self.constraintFile, self.cavityLigandFile, self.protein, n_processes = n_cores) #Dock in parallel
        
        #Compute distance to pharmacophoric point
        self.multiDistances = self.docking.assessAllDocks(self.multiDocks, self.hMulti, single = False)
    
    
        #######IDENTIFY QUASI-ACTIVES#########
        
        #Select up to 5 molecules where the pharmaocphore is within 2A of the pharmacophoric point for all specified pharmacophoric points
        self.multiQuasiActives = self.multiDistances.loc[self.multiDistances['distance'] < 2].drop_duplicates('smiles').head(5)
        
        
        #########REFINEMENT###############
        self.pharmElabsTestMulti = pd.DataFrame()

        with open(self.frag, 'r') as f:
            fragSmiles = f.read()

        
        self.multiPharmElabs = pharmElabs()
        self.multiPharmElabs.profileElabs = self.elaborate.makePharmElabsQuasiActives(self.multiQuasiActives, fragSmiles, filterMols = False) #make elaborations using the quasi actives profiles
        self.pharmElabsTestMulti = self.pharmElabsTestMulti.append(self.multiPharmElabs.profileElabs)

        #Now we have all of the molecules sampled from the quasi active profiles
        self.pharmElabsTestMultiSample = self.pharmElabsTestMulti.sample(n = numElabsRefinement, random_state = 10) #Sample the number of elaborations we want
        self.pharmElabsTestMultiFilter = self.elaborate.filterGeneratedMols(self.pharmElabsTestMultiSample, n_cores = n_cores)

        #Now prepare for docking 
        self.pharmElabsTestMultiFName = f'{self.storeLoc}/pharmsElabsTestMulti.sdf'
        self.pharmElabsMultiCountsIdx, self.pharmElabsMultiCountsSDFIdx = self.docking.getSDFs(self.pharmElabsTestMultiFilter, self.fragCore, self.pharmElabsTestMultiFName) 

        #Do docking
        self.pharmElabsTestMultiDocks, self.pharmElabsTestMultiFS = self.docking.dockLigandsMP(self.pharmElabsTestMultiFName, self.constraintFile, self.cavityLigandFile, self.protein, returnFitnessScore = True, n_processes = n_cores)

        #Compute ligand Efficiency
        self.pharmElabsTestMultiLigEff = self.docking.ligandEfficiency(self.pharmElabsTestMultiDocks, self.pharmElabsTestMultiFS)
        
        #Write the docks to file with the ligand efficiency as an attribute
        self.pharmElabsTestMultiDockedFName = f'{self.storeLoc}/pharmsElabsTestMultiDocked.sdf'
        w = Chem.SDWriter(self.pharmElabsTestMultiDockedFName)
        
        for idx, m in enumerate(self.pharmElabsTestMultiDocks):
            m.SetProp('STRIFE_LigEff_Score', str(self.pharmElabsTestMultiFS[idx]/m.GetNumHeavyAtoms()))
            w.write(m)

        w.close()
    
        #Standardise the name of the final df
        self.rankedElaborationsFinal = self.pharmElabsTestMultiLigEff

        if self.writeFinalElabs:
            self.rankedElaborationsFinal.to_csv(f'{self.storeLoc}/rankedElaborationsFinal.csv')

    
class HotspotSingle:
    def __init__(self, hotspotInfo, frag, clf, constraintFile, cavityLigandFile, protein):

        self.hotspots = hotspotInfo
        self.frag = frag
        self.clf = clf
        self.constraintFile = constraintFile
        self.cavityLigandFile = cavityLigandFile
        self.protein = protein

    
    def vector_angle(self, x, y):
        #Returns the angle between two numpy arrays
        cos_theta = np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))

        return np.arccos(cos_theta)

    def vector_distance(self, x, y):
        #Returns the distance between two numpy arrays
        diff = np.subtract(x,y)
        return np.linalg.norm(diff)


    def determineProfiles(self):
        #Essentially we just want to determine how many atoms we want to add,
        #in the case where we have an aromatic and the case where we don't

        #In the case where we don't have any aromatics, we just use the predicted atom count (which is the path - 1)

        #The case with aromatics is more difficult but we set the elab length as the combined distance plus 2-3 atoms
        #To allow for different placements of the aromatic ring

        #Fix donor and acceptor Counts
        acceptorCount = 0
        donorCount = 0


        if self.hotspots['type'] == 'Donor':
            donorCount += 1
        elif self.hotspots['type'] == 'Acceptor':
            acceptorCount += 1

        elabLengths = []
        profiles = []

        ###Profile Without Aromatic

        #Obtain number of atoms needed to get to the pharmacophoric point:
        dfForModelPrediction = pd.DataFrame({0:[0], 1:[0], 2:[1], 'dist': [self.hotspots['distFromExit']], 'ang': [self.hotspots['angFromExit']], 'aromaticEnRoute':[0]})
        predPathLength = [int(s) for s in self.clf.predict(dfForModelPrediction)][0]

        elabLengths.append(predPathLength - 1) #-1 because we predict path lengths, not the number of atoms required
        profiles.append([0, 0, acceptorCount, donorCount, 0])


        ###Profile With Aromatic
        dfForModelPredictionAro = pd.DataFrame({0:[0], 1:[0], 2:[1], 'dist': [self.hotspots['distFromExit']], 'ang': [self.hotspots['angFromExit']], 'aromaticEnRoute':[1]})
        predPathLengthAro = [int(s) for s in self.clf.predict(dfForModelPredictionAro)][0]

        elabLengths.append(predPathLengthAro + 1) #Predicted atom count +2
        elabLengths.append(predPathLengthAro + 2) #Predicted atom count +3
        elabLengths.append(predPathLengthAro + 3) #Predicted atom count +4

        profiles.append([0, 0, acceptorCount, donorCount, 1])
        profiles.append([0, 0, acceptorCount, donorCount, 1])
        profiles.append([0, 0, acceptorCount, donorCount, 1])

        return elabLengths, profiles
    
    
    
class HotspotMulti:

    def __init__(self, hotspotInfo, frag, clf, constraintFile, cavityLigandFile, protein):

        self.hotspots = hotspotInfo
        self.frag = frag
        self.clf = clf
        self.constraintFile = constraintFile
        self.cavityLigandFile = cavityLigandFile
        self.protein = protein

    def vector_angle(self, x, y):
        #Returns the angle between two numpy arrays
        cos_theta = np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))

        return np.arccos(cos_theta)

    def vector_distance(self, x, y):
        #Returns the distance between two numpy arrays
        diff = np.subtract(x,y)
        return np.linalg.norm(diff)


    def determineProfiles(self):
        #The Heuristic for determining the number of atoms we want is based on the
        #distance from the exit point to the furthest atom away, which has index 0
        
        #To account for side chains and aromatic groups, we specify a variety of different elaboration lengths:
            
            #For no aromatic - The predicted distance to the furthest pharmacophoric point plus 1 atom and 2 atoms
            #For aromatic - The predicted distance to the furthest pharmacophoric point plus 2-5 atoms



        #Set donor and acceptor Counts
        acceptorCount = 0
        donorCount = 0


        for k in self.hotspots.keys():
            if self.hotspots[k]['type'] == 'Donor':
                donorCount += 1
            elif self.hotspots[k]['type'] == 'Acceptor':
                acceptorCount += 1

        elabLengths = []
        profiles = []

        ###Profile Without Aromatic

        #Obtain number of atoms needed to get to the first hotspot and from there to the second hotspot:

        dfForModelPrediction = pd.DataFrame({0:[0], 1:[0], 2:[1], 'dist': [self.hotspots[0]['distFromExit']], 'ang': [self.hotspots[0]['angFromExit']], 'aromaticEnRoute':[0]})
        predPathLength = [int(s) for s in self.clf.predict(dfForModelPrediction)][0]

        elabLengths.append(predPathLength - 1) #-1 because we predict path lengths, not the number of atoms required
        elabLengths.append(predPathLength)
        elabLengths.append(predPathLength + 1)
        
        profiles.append([0, 0, acceptorCount, donorCount, 0])
        profiles.append([0, 0, acceptorCount, donorCount, 0])
        profiles.append([0, 0, acceptorCount, donorCount, 0])

        
        ###Profile With Aromatic
        dfForModelPredictionAro = pd.DataFrame({0:[0], 1:[0], 2:[1], 'dist': [self.hotspots[0]['distFromExit']], 'ang': [self.hotspots[0]['angFromExit']], 'aromaticEnRoute':[1]})
        predPathLengthAro = [int(s) for s in self.clf.predict(dfForModelPredictionAro)][0]

        elabLengths.append(predPathLengthAro + 1) #Predicted atom count +2
        elabLengths.append(predPathLengthAro + 2) #Predicted atom count +3
        elabLengths.append(predPathLengthAro + 3) #Predicted atom count + 4
        elabLengths.append(predPathLengthAro + 4) #Predicted atom count + 5


        profiles.append([0, 0, acceptorCount, donorCount, 1])
        profiles.append([0, 0, acceptorCount, donorCount, 1])
        profiles.append([0, 0, acceptorCount, donorCount, 1])
        profiles.append([0, 0, acceptorCount, donorCount, 1])


        return elabLengths, profiles

        


class pharmElabs:
    #Class to store things related to making elaborations on the quasi actives
    #We just initialise it as an empty class
    pass
    

if __name__=='__main__':
    
    parser = argparse.ArgumentParser()

    
    parser.add_argument('--fragment_sdf', '-f', type = str, required = True,
                        help = 'Location of fragment SDF. Can be an SDF file of a larger ligand for which the fragment is a substructure')
    parser.add_argument('--fragment_smiles', '-s', type = str, default = None,
                        help = 'Location of file containing fragment smiles string. The exit vector should be denoted by a dummy atom. Either you must specify an exit_vector_index or provide a fragment_smiles.')
    parser.add_argument('--exit_vector_idx', '-v', type = int, default = None,
                        help = 'The atomic index of the atom you wish to use as an exit vector. Either you must specify an exit_vector_index or provide a fragment_smiles. If you provide an exit_vector_index, then STRIFE can only use the molecule provided in fragment_sdf to make elaborations, whereas you can use fragment_smiles to make elaborations to substructures of the molecule provided in fragment_sdf.'  )
    
    
    parser.add_argument('--protein', '-p', type = str, required = True,
                        help = 'Location of protein pdb file (should have already been protonated)')
    parser.add_argument('--output_directory', '-o', type = str, default = '.', 
                        help = 'Directory to store output (default = current directory)')

    
    parser.add_argument('--hotspots_output', '-z', type = str, default = None,
                        help = 'Location of zip file containing hotspots output, if already calculated')
    parser.add_argument('--calculate_hotspots', '-c', type = str, default = None,
                        help = 'Location to save a FHM - if not None, then STRIFE will calculate an FHM, store it in the specified location and use it to generate molecules')
    #parser.add_argument('--calculated_hotspots_dir', '-d', type = str, default = None, 
    #                    help = 'Directory to store the calculated hotspots map - should only be used if calculate_hotspots = True')
    parser.add_argument('--load_specified_pharms', '-m', action = "store_true", 
                        help = 'Use pharmacophores that have been manually specfied instead of ones derived from FHMs. If True, the output_directory should contain at least one of donorHotspot.sdf or acceptorHotspot.sdf')    
    parser.add_argument('--path_length_model', type = str, default = 'models/pathLengthPred_saved.pickle', 
                        help = 'Location of saved SVM for predicting path distances')
    
    parser.add_argument('--model_type', '-t', type = int, default = 0,
                        help = 'Specify which setting you wish to generate the molecules with: 0 -> default STRIFE algorithm, 1 -> simultaneously satisfy multiple pharmacophoric points (only recommended if you have manually specified the pharmacophores), 2 -> run the STRIFE algorithm without refinement. Default %(default)s')
    parser.add_argument('--number_elaborations', '-n', type = int, default = 250,
            help = 'Final number of elaborations for the model to generate. Default: %(default)s')
    parser.add_argument('--number_elaborations_exploration', '-e', type = int, default = 250,
            help = 'Number of elaborations to make per pharmacophoric point in the exploration phase. Default: %(default)s')
    parser.add_argument('--name', type = str, default = None,
                        help = 'Model name for saving the model. If None (the default argument) then the model will be saved as STRIFE_{date and time}.pickle')
    parser.add_argument('--write_elaborations_dataset', '-w', action = "store_true", 
            help='Save the DataFrame containing the final elaborations generated by STRIFE as rankedElaborationsFinal.csv')

    parser.add_argument('--num_cpu_cores', '-cpu', type = int, default = 1, 
            help='Number of CPU cores to use for docking and other computations. Specifiying -1 will use all available cores')

    #TODO
    #parser.add_argument('--compute_hotspot_distance', action = "store_true",
            #help='Optional flag which will compute the distance of ligand pharmacophores to the nearest pharmacophoric point and save as a molecule property')

    arguments = parser.parse_args()
    
    #Define STRIFE model
    STRIFE_model = STRIFE(arguments)
    
    #Set it running
    STRIFE_model.run(arguments)
    
    #Pickle the STRIFE Class and SAVE
    if arguments.name is None:
        run_id = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        
        with open(f'{arguments.output_directory}/STRIFE_{run_id}.pickle', 'wb') as f:
            pickle.dump(STRIFE_model, f)
    else:
        with open(f'{arguments.output_directory}/STRIFE_{arguments.name}.pickle', 'wb') as f:
            pickle.dump(STRIFE_model, f)
    

    
    
