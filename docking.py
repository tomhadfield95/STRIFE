#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 16:22:22 2021

@author: Tom Hadfield

Class for docking elaborations using GOLD and doing operations with docked structures

"""

#import elabHotspotsBase
#import elabHotspotsCCDC

#########Standard Libraries##########
import json
import numpy as np
from numpy import linalg
from scipy import stats
import pandas as pd
import multiprocessing as mp #For parallelising
import time
import matplotlib.pyplot as plt
import sys
from random import sample
from functools import partial
import glob

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
from rdkit.Chem import rdchem
from rdkit import DataStructs
from rdkit.Chem.Draw import IPythonConsole
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') #Supress annoying RDKit output
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures


import pickle
import addPharmProfile
import os

from ccdc.docking import Docker
from ccdc.io import MoleculeReader, EntryReader
import tempfile


import sys
from argparse import ArgumentParser
from platform import platform
from pathlib import Path
from os import mkdir, chdir
from shutil import rmtree, copy as cp
from time import time
from dataclasses import dataclass, field
from multiprocessing import Pool


#Create feature factory
fdefName = os.path.join(RDConfig.RDDataDir,'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

class docking:
    
    def __init__(self):
        pass
    
    
    def ligandEfficiency(self, mols, scores):
        #Calculate the ligand efficiency as docking_score/number_heavy_atoms

        smiles = []
        le = []

        for i, m in enumerate(mols):
            smiles.append(Chem.MolToSmiles(m))
            le.append(scores[i]/m.GetNumHeavyAtoms()) 

        leDF = pd.DataFrame({'smiles':smiles, 'ligEfficiency': le}).sort_values('ligEfficiency', ascending = False)
        return leDF
    
    def vectorDistance(self, p1, p2):
        return np.linalg.norm(p1 - p2)
    
    def assessAllDocks(self, mols, hotspot, single = True):
        smiles = []
        dist = []
        for m in mols:
            smiles.append(Chem.MolToSmiles(m))
            dist.append(self.assessSingleDock(m, hotspot, single))
        
        distanceDF = pd.DataFrame({'smiles':smiles, 'distance':dist}).sort_values('distance', ascending = True)
        return distanceDF
    
    
    
    def assessSingleDock(self, mol, hotspot, single = True):
        #Compute the distance between the pharmacophores in a mol and the corresponding hotspot
        

        if single == True:
            #i.e. we're comparing to a single hotspot
            #Check how close the pharmacophores in the docked molecule are to the hotspot of interest
            distanceToPharm = []
            feats = factory.GetFeaturesForMol(mol)
            for feat in feats:
                if feat.GetFamily() == hotspot['type']:
                    #Compute Distance to the  hotspot and take the pharmacophore with the smallest distance
                    pharmPosition = np.array(mol.GetConformer().GetAtomPosition(feat.GetAtomIds()[0]))
                    distanceToPharm.append(self.vectorDistance(pharmPosition, hotspot['position']))
            return min(distanceToPharm)
        else:
            #i.e. we're assessing closeness to multiple hotspots
            
            distances = []
            feats = factory.GetFeaturesForMol(mol)

            for k in hotspot.keys():
                pharmPoint = hotspot[k]
                distanceToPharm = []

                for feat in feats:
                    if feat.GetFamily() == pharmPoint['type']:
                        #Compute Distance to the hotspot and take the pharmacophore with the smallest distance
                        pharmPosition = np.array(mol.GetConformer().GetAtomPosition(feat.GetAtomIds()[0]))
                        distanceToPharm.append(self.vectorDistance(pharmPosition, pharmPoint['position']))
                if len(distanceToPharm) == 0:
                    #i.e. no matching pharmacophores
                    distances.append(100) #Append really big distance
                elif len(distanceToPharm) > 0:
                    distances.append(min(distanceToPharm))

            return max(distances) #Return the maximum of the distances to the hotspot centres


            '''
            distances = []
            for i in [1,2]:
                hotspotOfInterest = hotspot[i]
                feats = factory.GetFeaturesForMol(mol)
                distanceToPharm = []
                for feat in feats:
                   if feat.GetFamily() == hotspotOfInterest['type']:
                        #Compute Distance to the three hotspot and take the pharmacophore with the smallest distance
                        pharmPosition = np.array(mol.GetConformer().GetAtomPosition(feat.GetAtomIds()[0]))
                        distanceToPharm.append(self.vector_distance(pharmPosition, hotspotOfInterest['position']))

            if len(distanceToPharm) == 0:
                #i.e. there were no matching pharmacophores
                distances.append(100) #Append really big distance
            elif len(distanceToPharm) > 0:
                distances.append(min(distanceToPharm))
            return max(distances) #Return the maximum of the distances to the hotspot centres
            '''
    
    
    
    
    
    def flexDocking(self, ligandsToDockSDF, constraintFile, cavityLigandFile, protein, n_docks = 10, outputFile = 'docked_ligands.sdf', returnFitnessScore = False, conf_file = None, trackFlexibility = None):

        #for trackFlexibility - if not None, then should be a list: [resID, chiNumber] - #e.g. ['TYR119', 'chi2']
            
            
        if conf_file is None:
            return 'Please specify a conf file that supports flexible docking'

        else:
            #i.e. we use a conf_file to set up the docking (if we're going to be doing something like flexible docking

            #First we need to change the conf_file so that it has the right ligand file so we dock the correct ligands

            with open(conf_file, 'r') as f:
                cf = f.readlines()

            cf_new = []

            for l in cf:
                if 'ligand_data_file' in l:
                    cf_new.append(f'ligand_data_file {ligandsToDockSDF} {n_docks}\n')
                else:
                    cf_new.append(l)

            with open(conf_file, 'w') as f:
                for l in cf_new:
                    f.write(l)

            #Now the conf_file should have the correct ligand information but everything else should be the same

            #Set up docking with conf file
            settings = Docker.Settings.from_file(conf_file)
            docker = Docker(settings = settings)

        #Set scaffold constraint:
        scaffold =  MoleculeReader(constraintFile)[0]
        settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))


        settings.fitness_function = 'plp'
        settings.autoscale = 10.
        settings.early_termination = False
        batch_tempd = tempfile.mkdtemp()
        settings.output_directory = batch_tempd
        settings.output_file = outputFile


        '''
        #Define protein, binding site, scaffold and ligands to dock
        protein = settings.proteins[0]
        crystal_ligand = MoleculeReader(cavityLigandFile)[0]
        crystal_ligand.identifier = 'crystalLigand' #Change the identifier of the cavity ligand
        settings.binding_site = settings.BindingSiteFromLigand(protein, crystal_ligand, 10.0) #Define the binding site

        ligandsToDock = MoleculeReader(ligandsToDockSDF)
        settings.add_ligand_file(ligandsToDockSDF, 10) #Generate 10 poses per 

        #Set scaffold constraint:
        scaffold =  MoleculeReader(constraintFile)[0]
        settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))
        '''

        #Dock the ligands
        results = docker.dock()

        #Get docked poses:
        docksLocation = results.ligands.file_name[0]
        topRanked = glob.glob(f'{docksLocation[0:docksLocation.index(outputFile)]}ranked*_m*_1.sdf')
        mols = [Chem.SDMolSupplier(fn)[0] for fn in topRanked]
        
        
        if trackFlexibility is not None:
            
            greatest_torsion_changes = []
            
            
            #Get Mol Indices
            for i in range(len(topRanked)):
                docksForOneLigand = glob.glob(f'{docksLocation[0:docksLocation.index(outputFile)]}ranked*_m{i + 1}_*.sdf')
                torsion_changes = []
                
                for dock in docksForOneLigand:
                    rt = EntryReader(dock)[0].attributes['Gold.Protein.RotatedTorsions']
                    
                    for l in rt.split('\n'):
                        if trackFlexibility[0] in l and trackFlexibility[1] in l:
                            end_torsion = float(l[l.index('final') + 5:l.index('input')])
                            start_torsion = float(l[l.index('input') + 5:l.index('|library')])
                    torsion_changes.append(np.abs(end_torsion - start_torsion))
                
                greatest_torsion_changes.append(max(torsion_changes))
        
        
        locs = glob.glob(f'{docksLocation[0:docksLocation.index(outputFile)]}ranked*_m*.sdf')

        #Either return just the docks or the docks and associated fitness scores

        if returnFitnessScore and trackFlexibility is None:
            scores = [float(EntryReader(fn)[0].attributes['Gold.PLP.Fitness']) for fn in topRanked] #Get PLP scores for each of the top-ranked docks
            return mols, scores
        elif returnFitnessScore and trackFlexibility is not None:
            scores = [float(EntryReader(fn)[0].attributes['Gold.PLP.Fitness']) for fn in topRanked] #Get PLP scores for each of the top-ranked docks
            return mols, scores, greatest_torsion_changes
        else:
            return mols



    def dockLigands(self, ligandsToDockSDF, constraintFile, cavityLigandFile, protein, n_docks = 10, outputFile = 'docked_ligands.sdf', returnFitnessScore = False):
        
        #Method for the (rigid) docking of ligands, without Multiprocessing
        #For faster results, it is recommended that you use dockLigandsMP()
        #To perform docking with a flexible side chain, we recommend that you use flexDocking()
        
        
        
        #Set up docking class
        docker = Docker()
        settings = docker.settings
        settings.add_protein_file(protein)
        protein = settings.proteins[0]

        crystal_ligand = MoleculeReader(cavityLigandFile)[0]
        crystal_ligand.identifier = 'crystalLigand' #Change the identifier of the cavity ligand
        settings.binding_site = settings.BindingSiteFromLigand(protein, crystal_ligand, 10.0) #Define the binding site

        ligandsToDock = MoleculeReader(ligandsToDockSDF)
        settings.add_ligand_file(ligandsToDockSDF, 10) #Generate 10 poses per ligand

        #Set scaffold constraint:
        scaffold =  MoleculeReader(constraintFile)[0]
        settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))

        settings.fitness_function = 'plp'
        settings.autoscale = 10.
        settings.early_termination = False
        batch_tempd = tempfile.mkdtemp()
        settings.output_directory = batch_tempd
        settings.output_file = outputFile

        #Dock the ligands
        results = docker.dock()

        #Get docked poses:
        docksLocation = results.ligands.file_name[0]
        topRanked = glob.glob(f'{docksLocation[0:docksLocation.index(outputFile)]}ranked*_m*_1.sdf')
        mols = [Chem.SDMolSupplier(fn)[0] for fn in topRanked]
        
        #Either return just the docks or the docks and associated fitness scores

        if returnFitnessScore:
            scores = [float(EntryReader(fn)[0].attributes['Gold.PLP.Fitness']) for fn in topRanked] #Get PLP scores for each of the top-ranked docks
            return mols, scores
        else:
            return mols


    def dockLigandsMP(self, ligandsToDockSDF, constraintFile, cavityLigandFile, protein, n_processes = 7, ndocks = 10, outputFile = 'docked_ligands.sdf', returnFitnessScore = False):

        """
        Dock the molecules from the supplied input file in parallel.
        Adapted from the gold_multi_map.py file supplied by ccdc
        """

        t0 = time()  # Script start time


        if not n_processes > 0:
            print(f"Error! Number of processes must be an integer greater than zero.")
            sys.exit(1)

        #Get number of molecules to be docked
        with EntryReader(ligandsToDockSDF) as reader:
            n_molecules = len(reader)

        print(f"There are {n_molecules} molecules to dock on {n_processes} processes...")

        output_dir = tempfile.mkdtemp() #Create temporary output file

        #########################################

        # Determine the sets of parameters defining the chunks...
        chunks = []  # List of records that define the chunks

        # Work out the size of each chunk, and hence the start and finish indices in the input file...
        basic_size = n_molecules // n_processes  # Basic size of a chunk, which must obviously be integral
        remainder  = n_molecules %  n_processes  # Number of molecules that would not be included in basic-sized chunks

        finish = 0  # Finish index

        for n in range(1, n_processes + 1):  # Recall that the number of chunks is the same as the number of processes
            start = finish + 1  # Start index
            chunk_size = basic_size + 1 if n <= remainder else basic_size  # Add one to the basic chunk sizes until the remainder are taken care of
            finish = start + chunk_size - 1
            chunks.append(chunkForMP(n=n, start=start, finish=finish, output_dir=output_dir, crystal_ligand = cavityLigandFile, constraint = constraintFile, protein = protein, ligand_file = ligandsToDockSDF, ndocks=ndocks))
            print(f"chunk {n}: size: {chunk_size}; start: {start}, finish: {finish}.")

        ##########################################

        # Dock the chunks in parallel...
        with Pool(n_processes) as pool:
            _ = pool.map(do_chunk, chunks)  # No output; docks are saved in the output_dir directory

        #print(f"Finished docking in {time() - t0:.1f} seconds.")

        #Now return the top ranked mols
        topRanked = sorted(glob.glob(f'{output_dir}/*/ranked*_m*_1.sdf')) #sort the list for consistency
        mols = [Chem.SDMolSupplier(fname)[0] for fname in topRanked]
        
        # All done.
        print(f"Finished in {time() - t0:.1f} seconds.")
        
        if returnFitnessScore:
            scores = [float(EntryReader(fn)[0].attributes['Gold.PLP.Fitness']) for fn in topRanked] #Get PLP scores for each of the top-ranked docks
            return mols, scores
        else:
            return mols


    def prepareForDocking(self, Hotspot, core, fname):
        
        #Hotspot is a HotspotSingle or HotspotMulti

        Hotspot.dockingFname = fname #Save the name of the file we'll save the hotspots to. 
        Hotspot.countsIdx, Hotspot.countsSDFIdx = self.getSDFs(Hotspot.profileElabs, core, fname) #Saves sdfs to fname
        
        return Hotspot


    def getMultipleSDFs(self, data, core, outNameCore):
    
        #Iterate over the different indices (i.e. the different pharmacophoric profiles) and pass to getSDFs
        
        for idx in data['idx'].drop_duplicates():
            d = data.loc[data['idx'] == idx] 
            self.getSDFs(d, core, f'{outNameCore}_p{idx}.sdf')
        
    
    def getSDFs(self, data, core, outName):
        #Core is the embedded mol we use for the constrained embedding
        
        mols = []
        dfIdx = []
        sdfIdx = []
        sdfIdxNum = 0
        for idx, g in enumerate(list(data['gen'])):
            dfIdx.append(idx)
            m = Chem.MolFromSmiles(g)
            
            try:
                AllChem.ConstrainedEmbed(m, core)
                fail = 0
            except:
                fail = 1
                print(idx)
    
    
            #Now addHs and reembed:
            mHs = Chem.AddHs(m)
            try:
                AllChem.ConstrainedEmbed(mHs, m)
                fail = 0
            except:
                fail = 1
                print(idx)
            
            if fail == 0:
                mols.append(mHs)
            else:
                mols.append('Could not embed molecule - not docking')
        

        #Now write mols to file
        w = Chem.SDWriter(outName)
        for m in mols:
            if m != 'Could not embed molecule - not docking':
                w.write(m)
                sdfIdx.append(sdfIdxNum)
                sdfIdxNum += 1
            else:
                sdfIdx.append(-1) #Indicates that we weren't able to dock this molecule
            
        return dfIdx, sdfIdx


####Additional Code for parallel docking####

@dataclass
class chunkForMP:

    #Contains the information to provide to GOLD for docking
    #Store a bunch of chunks in a list and use them to do the docking in parallel

    n: int       # Chunk number
    start: int   # Index of first molecule in chunk
    finish: int  # Index of last molecule in chunk
    ndocks: int  # Number of docking poses to generate
    #conf_file: Path   # GOLD configuration file
    output_dir: Path  # Output dir, in which the chunk sub-directory will be created
    crystal_ligand: Path #Location of the crystal ligand (used for binding site definition)
    constraint: Path #Location of mol2 file used to constrain the docking
    protein: Path #Location of protein pdb File
    ligand_file: Path #Location of ligands to be docked

    dir: Path = field(init=False) # Sub-directory for chunk, see __post_init__ below

    def __post_init__(self):
        self.dir = f'{self.output_dir}/chunk_{self.n:02d}'


def do_chunk(chunk):
    
    """
    Dock a chunk of the input file.

    :param chunk: a record holding the parameters defining the chunk

    As we can't return a GOLD results object from a pool process (it can't be pickled as it wraps C++ objects),
    we simply return a boolean recording whether GOLD exited with a 'success' status code.
    """

    
    #Set up docking class
    docker = Docker()
    settings = docker.settings
    settings.add_protein_file(chunk.protein)
    
    settings.fitness_function = 'plp'
    settings.autoscale = 10.
    settings.early_termination = False
#     batch_tempd = tempfile.mkdtemp()
#     settings.output_directory = batch_tempd
#     settings.output_file = outputFile
    
    
    # Create and enter the sub-directory for this chunk...
    mkdir(chunk.dir)
    chdir(chunk.dir)
    settings.output_directory = '.'  # Ensure GOLD writes output to the chunk sub-directory

    
    #Define protein, binding site, scaffold and ligands to dock
    protein = settings.proteins[0]
    crystal_ligand = MoleculeReader(chunk.crystal_ligand)[0]
    crystal_ligand.identifier = 'crystalLigand' #Change the identifier of the cavity ligand
    settings.binding_site = settings.BindingSiteFromLigand(protein, crystal_ligand, 10.0) #Define the binding site
    
    # Specify the chunk of molecules to dock...
#     ligand_file = settings.ligand_files[0]  # The ligand file info will be overwritten, so store for reference below
    settings.clear_ligand_files() #Clear just to be sure
    settings.add_ligand_file(chunk.ligand_file, ndocks=chunk.ndocks, start=chunk.start, finish=chunk.finish)

    #Set scaffold constraint:
    scaffold =  MoleculeReader(chunk.constraint)[0]
    settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))
    
    
    # Run docking...
    print(f"Starting chunk {chunk.n} (ligand indices {chunk.start} - {chunk.finish})...")
    docker = Docker(settings=settings)
    results = docker.dock()

    
    
    print(f"Finished chunk {chunk.n}.")
    # As we can't return the results (as they are not picklable) and the poses have already been written to disk, we just return the status code

    return results.ligands.file_name[0] #return docks location
