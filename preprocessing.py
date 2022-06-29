#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:12:51 2021

@author: hadfield
"""

# Generic imports
from pathlib import Path
from shutil import make_archive, rmtree, unpack_archive
import numpy as np
import pandas as pd
import sys
from itertools import product


#rdkit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMMPA


#CCDC imports
from ccdc.protein import Protein

# Hotspots imports
from hotspots.atomic_hotspot_calculation import _AtomicHotspotResult
from hotspots.calculation import Runner
from hotspots.grid_extension import Grid
from hotspots import hs_io


class preprocessing:
    
    #Class for various preprocessing tasks associated with the STRIFE algorithm.
    
    def __init__(self):
        pass
    
    ###Hotspot Calculation###
    #Run calculateFHM() to calculate and store a Fragment Hotspot Map
    def calculateFHM(self, pdbFile, outputDirectory = ''):
        #pdbfile - Location of a pdb file that you want to calculate the hotspots for
        #outputDirectory - Directory where the hotspots output will be stored. (Default is to store in the current working directory)
    
        prot_path = Path(Path.cwd(), pdbFile)
        out_path = Path(Path.cwd(), outputDirectory)
        if not out_path.exists(): out_path.mkdir()
    
        self.run_hotspot_calculation(protein_path=prot_path,
                                out_dir=out_path,
                                method="ghecom")

    
    
    def prepare_protein(self, protein_path):
        prot = Protein.from_file(str(protein_path))
        prot.remove_all_waters()
        prot.remove_all_metals()
        prot.add_hydrogens()
        for l in prot.ligands:
            prot.remove_ligand(l.identifier)
        return prot

    def save_superstar_grids(self, out_path, hs_runner):
        """
        Saves and Zips the SuperStar grids from the hotspots.calculation.Runner.
        Need to also save the buriedness grid.
        :param hs_runner:
        :return:
        """
        superstar_dir = Path(out_path, "superstar_grids")
        if not superstar_dir.exists(): superstar_dir.mkdir()
    
        # Need to save the buriedness grid as well as the superstar grids
        hs_runner.buriedness.write(str(Path(superstar_dir, "buriedness.ccp4").resolve()))
        for s in hs_runner.superstar_grids:
            s.grid.write(str(Path(superstar_dir, "superstar_{}.ccp4".format(s.identifier)).resolve()))
        make_archive(superstar_dir, 'zip', superstar_dir)
        rmtree(superstar_dir)
    
    def create_atomic_hotspots(self, superstar_grids_dir):
        """
    
        :param superstar_grids_dir: path to where the superstar grids are stored
        :return:
        """
        atomic_hotspots = []
    
        # Read in the SuperStar and Buriedness info
        probes = ['donor', 'acceptor', 'apolar', 'positive', 'negative']
        b_grid = Grid.from_file(str(Path(superstar_grids_dir, 'buriedness.ccp4').resolve()))
    
        for p in probes:
            g_path = Path(superstar_grids_dir, f'superstar_{p}.ccp4')
            if g_path.exists():
                print(f" found grid for probe of type {p}")
                p_grid = Grid.from_file(str(g_path.resolve()))
                ahs = _AtomicHotspotResult(identifier=p, grid=p_grid, buriedness=b_grid)
                atomic_hotspots.append(ahs)
            else:
                continue
    
        return atomic_hotspots
    
    def run_hotspot_calculation(self, protein_path, out_dir, method="ghecom"):
        """
        Runs the hotspots calculation on the specified PDB structure
        :return:
        """
        h = Runner()
        settings = h.Settings()
        settings.nrotations = 3000
        settings.apolar_translation_threshold = 15
        settings.polar_translation_threshold = 15
        settings.sphere_maps = False
    
    
        # Check if SuperStar and Ghecom have already been run.
        super_archive_path = Path(out_dir.parent, "superstar_grids.zip")
    
        if super_archive_path.exists():
            super_tmp_path = Path(out_dir.parent, super_archive_path.stem)
    
            if not super_tmp_path.exists(): super_tmp_path.mkdir()
            unpack_archive(super_archive_path, super_tmp_path, 'zip')
            b_grid = Grid.from_file(str(Path(super_tmp_path, 'buriedness.ccp4').resolve()))
    
            result = h.from_superstar(protein=self.prepare_protein(protein_path),
                                            superstar_grids=self.create_atomic_hotspots(super_tmp_path),
                                            buriedness=b_grid,
                                            charged_probes=False,
                                            settings=settings,
                                            clear_tmp=True)
            rmtree(super_tmp_path)
    
        else:
    
            result = h.from_protein(protein=self.prepare_protein(protein_path),
                                    charged_probes=False,
                                    probe_size=7,
                                    buriedness_method=method,
                                    cavities=None,
                                    nprocesses=3,
                                    settings=settings)
    
            # Save and zip the SuperStar Grids:
            self.save_superstar_grids(out_dir, h)
    
        # Save and zip the Results
        #with hs_io.HotspotWriter(str(out_dir), visualisation="pymol", grid_extension=".ccp4", zip_results=True) as writer:
        with hs_io.HotspotWriter(str(out_dir), grid_extension=".ccp4", zip_results=True) as writer:
            writer.write(result)
    
        print(f"out_file: {str(Path(out_dir, 'out.zip'))}")
    
        return Path(out_dir, 'out.zip')
    

    
    ####Fragment Preprocessing####
    
    def preprocessFragment(self, fragSmiles, mol3D):
    
        #Take a smiles (with a dummy atom denoting the exit vector) and do the processing required to make elaborations
        #1. Create a molecule of fragments with the 3d coordinates extracted from the sdfFile (which can be a molecule of the fragment or some larger ligand)
        #2. Align the fragment molecule to the molecule with the dummy atom, so that we can extract information about the exit vector and set up the constraint file for the docking
        
        molWithBondOrders = self.assignBondOrderTo3DMol(mol3D)
        
        
        fragMol = Chem.MolFromSmiles(fragSmiles)
        
        if molWithBondOrders.GetNumHeavyAtoms() > fragMol.GetNumHeavyAtoms():
            #The conformer we're using is a larger molecule than the fragment
            #either because we're paring back the fragment a little or because we're running on a test set
            #with 'real answers'
            
            #Compute ref elaboration
            refElab = self.returnGenElaborationWithDummy(fragSmiles, Chem.MolToSmiles(molWithBondOrders))
            
            #Get a 3D conformer of the frag
            fragConf = self.getFragConf(molWithBondOrders, refElab, fragSmiles)
            
        elif molWithBondOrders.GetNumHeavyAtoms() == fragMol.GetNumHeavyAtoms():
            
            fragConf = molWithBondOrders #make naming consistent with the other case
            
        
        #Now align fragMol (which is a 2D molecule with a dummy atom) with fragConf, so that the atom indices
        #are consistent
        
        #As a result of this we'll be able to determine which atom in fragConf is the exit vector
        #fragConf should be a substructure of fragMol (identical apart from the dummy atom)
        
        if fragMol.HasSubstructMatch(fragConf):
            
            #Need to reorder the atom indices of fragMol so that they mirror those of fragConf (with the dummy atom as the final idx)
            for atom in fragMol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    dummyIdx = atom.GetIdx()
    
            atomIndices = list(fragMol.GetSubstructMatch(fragConf))
            atomIndices.append(dummyIdx)
            
            fragMol = Chem.RenumberAtoms(fragMol, atomIndices)
            
            
            #Now we want to get the idx and 3D coordinates of the neighbour of the dummy atom
            
            for atom in fragMol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    exitVectorIdx = [nei.GetIdx() for nei in atom.GetNeighbors()][0] #Dummy atom only has one neighbour
            
            exitVectorPos = np.array(fragConf.GetConformer().GetAtomPosition(exitVectorIdx))
            
            
            
            
            #Finally we want to create a molecule we can use for the constraint in gold docking
            #Need to add Hydrogens but remove the ones attached to the exit vector
            
            fragConf2D = Chem.MolFromSmiles(Chem.MolToSmiles(fragConf))
            
            #Renumber atoms so they're the same as fragConf
            fragConf2D = Chem.RenumberAtoms(fragConf2D, fragConf2D.GetSubstructMatch(fragConf))
            
            fragConf2D = Chem.AddHs(fragConf2D)
            AllChem.ConstrainedEmbed(fragConf2D, fragConf) #Create 3D molecule
            
            fragConf2DRW = Chem.RWMol(fragConf2D)
            toDelete = []
            
            for n in fragConf2DRW.GetAtomWithIdx(exitVectorIdx).GetNeighbors():
                if n.GetSymbol() == 'H':
                    toDelete.append(n.GetIdx())
            
            for idx in sorted(toDelete, reverse = True):
                fragConf2DRW.RemoveAtom(idx)
            
            Chem.SanitizeMol(fragConf2DRW)
            
            return fragConf, exitVectorIdx, exitVectorPos, fragConf2DRW
            
        else:
            return -1 #Error
    
    
    def returnGenElaborationWithDummy(self, frag, gen):
        #Both inputs are smiles strings - frag should only have a single dummy atom
        #Return the elaboration
        #i.e. the part of gen which is not part of frag
    
        generated_elaboration = self.single_elaboration_2(frag, gen)
    
        if generated_elaboration == 'Unable to determine elaboration':
            generated_elaboration = self.single_elaboration(frag, gen)
    
        return generated_elaboration


    def single_elaboration_2(self, frag, full):
        #take a molecule and fragment of that molecule - return the elaborated structure by enumerating cuts
        #to the molecule and working with the pair that matches the fragment
    
        output_smiles = None
    
        cuts = self.fragment_mol_for_elaboration(full, full)
    
        #Now iterate through:
    
        for _, _, frag1, frag2 in cuts:
            frag1, frag2
            c1 = self.compare_by_inchi(frag1, frag)
            c2 = self.compare_by_inchi(frag2, frag)
    
            if c1 == 1:
                output_smiles = frag2
            if c2 == 1:
                output_smiles = frag1
    
        if output_smiles is None:
            return 'Unable to determine elaboration'
    
    
        return output_smiles
    
    def single_elaboration(self, frag, full):
        #Aligns by fragments, finds the exit node and fragments on that.
        #Returns the molecule which doesn't have a substruct match with frag
    
    
        #FRAG SHOULD ONLY HAVE A SINGLE DUMMY ATOM
        if len(frag.split('.')) > 1:
            frag = frag.split('.')[0]
    
    
        length_of_elab = Chem.MolFromSmiles(full).GetNumHeavyAtoms() - Chem.MolFromSmiles(frag).GetNumHeavyAtoms()
    
        if length_of_elab < 2:
            return '*'
        else:
    
            (aligned_mol, aligned_frag), to_keep, exit = self.align_smiles_by_frags(full, frag)
    
            #Get neighbours of exit atom
            neighbors = aligned_mol.GetAtomWithIdx(exit[0]).GetNeighbors()
    
            #Iterate through neighbours and take the one which is in the elaborated section
            for neighbor in neighbors:
                if neighbor.GetIdx() not in to_keep:
                    bond = aligned_mol.GetBondBetweenAtoms(aligned_mol.GetAtomWithIdx(exit[0]).GetIdx(), neighbor.GetIdx())
    
            new_mol = Chem.FragmentOnBonds(aligned_mol, [bond.GetIdx()])
    
            #We want to make sure we take the elaborated bit
            # Include dummy in query
            qp = Chem.AdjustQueryParameters()
            qp.makeDummiesQueries=True
            qfrag = Chem.AdjustQueryProperties(Chem.MolFromSmiles(frag),qp)
    
    
    
            #new_mol.HasSubstructMatch(qfrag)
            new_mol_smiles = Chem.MolToSmiles(new_mol).split('.')
            #return new_mol_smiles
    
            #Check the two fragments
            for smiles in new_mol_smiles:
                if Chem.MolFromSmiles(smiles).HasSubstructMatch(qfrag) == False:
                    return smiles
                else:
                    print(smiles)
    
    
            for smiles in new_mol_smiles:
                #Get rid of dummy atoms
                smiles2 = self.remove_dummy_atom(smiles)
                #Compute inchikey
                inchisplit = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smiles2))
    
    
                #Now get the inchi key of the original fragment
                frag2 = self.remove_dummy_atom(frag)
                inchifrag = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(frag2))
    
                if inchisplit != inchifrag:
                    return smiles


    def align_smiles_by_frags(self, smiles_mol, smiles_frag):
        #Amended function which takes a single fragment as input
        try:
            smiles_frags = smiles_frag + '.[*:2]'
            mols_to_align = [Chem.MolFromSmiles(smiles_mol), Chem.MolFromSmiles(smiles_frags)]
            frags = [Chem.MolFromSmiles(smiles_frag)]
    
            # Include dummy in query
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


    def fragment_mol_for_elaboration(self, smi, cid, pattern="[#6+0;!$(*=,#[!#6])]!@!=!#[*]"):
        mol = Chem.MolFromSmiles(smi)
    
        #different cuts can give the same fragments
        #to use outlines to remove them
        outlines = set()
    
        if (mol == None):
            sys.stderr.write("Can't generate mol for: %s\n" % (smi))
        else:
            frags = rdMMPA.FragmentMol(mol, minCuts=1, maxCuts=1, maxCutBonds=100, pattern=pattern, resultsAsMols=False)
            for _, chains in frags:
                output = '%s,%s,%s' % (smi, cid, chains)
                if (not (output in outlines)):
                    outlines.add(output)
            if not outlines:
                # for molecules with no cuts, output the parent molecule itself
                outlines.add('%s,%s,,' % (smi,cid))
    
    
            out_list = list(outlines)
            formatted = []
    
    
        for l in out_list:
            formatted.append(l.replace('.', ',').split(','))
    
        return formatted


    def compare_by_inchi(self, mol1, mol2):
        #take two smiles strings and assess if they equal by converting to inchi strings
        #helper function for single_elaboration_2
        #Remove dummy atoms
        mol1 = mol1.split('.')[0]
        mol2 = mol2.split('.')[0]
    
        mol1 = self.remove_dummy_atom(mol1)
        mol2 = self.remove_dummy_atom(mol2)
    
    
    
        inchi1 = Chem.MolToInchi(Chem.MolFromSmiles(mol1))
        inchi2 = Chem.MolToInchi(Chem.MolFromSmiles(mol2))
    
    
        if inchi1 == inchi2:
            return 1
        else:
            return 0

    def remove_dummy_atom(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol2 = AllChem.ReplaceSubstructs(mol, Chem.MolFromSmiles('*'), Chem.MolFromSmiles('[H]'), True)[0]
        mol3 = Chem.RemoveHs(mol2)
    
        return Chem.MolToSmiles(mol3)


    
    
    def assignBondOrderTo3DMol(self, m):
        #RDKit is behaving strangely - it's loading the pdbbind molecules in without bond orders
        #but still retains that information somehow 
        #This function takes in the molecule and returns the 3D molecule with the correct bond orders
        
        m2 = Chem.MolFromSmiles(Chem.MolToSmiles(m))
        
        if m.HasSubstructMatch(m2):
            #Renumber indices
            m_renum = Chem.RenumberAtoms(m, m.GetSubstructMatch(m2))
            conf = Chem.Conformer(m.GetNumHeavyAtoms())
            
            for i in range(m.GetNumHeavyAtoms()):
                pos = list(m_renum.GetConformer().GetAtomPosition(i))
                conf.SetAtomPosition(i, pos)
            m2.AddConformer(conf, assignId = True)
            return m2
        else:
            return 'Unable to find conformer in SDF File'
    
    
    
    def getFragConf(self, mol, smi_linker, smi_frag):
        """ aligns the 3D coordinates of the fragment to the
        coordinates of a reference molecule that contains the fragment.
        Arguments
          - mol: the reference molecule
          - smi_linker: the elaboration with exit vector
          - smi_frag: the fragments with exit vector
        """
    
        try:
                   
            smi_frags = smi_frag + '.[*:2]'
            ####################SHOULD ALL BE THE SAME#####################
            #frags = [Chem.MolFromSmiles(frag) for frag in smi_frags.split(".")]
            mol = Chem.RemoveHs(mol)
            frags = Chem.RemoveHs(Chem.MolFromSmiles(smi_frags))
            linker = Chem.RemoveHs(Chem.MolFromSmiles(smi_linker))
            # Include dummy in query
            qp = Chem.AdjustQueryParameters()
            qp.makeDummiesQueries=True
            # Renumber based on frags (incl. dummy atoms)
            aligned_mols = []
    
            sub_idx = []
            # Align to frags and linker
            qfrag = Chem.AdjustQueryProperties(frags,qp)
            frags_matches = list(mol.GetSubstructMatches(qfrag, uniquify=False))
            qlinker = Chem.AdjustQueryProperties(linker,qp)
            linker_matches = list(mol.GetSubstructMatches(qlinker, uniquify=False))
    
            ################################################################
    
            # Loop over matches
            for frag_match, linker_match in product(frags_matches, linker_matches):
                # Check if match
                f_match = [idx for num, idx in enumerate(frag_match) if frags.GetAtomWithIdx(num).GetAtomicNum() != 0]
                l_match = [idx for num, idx in enumerate(linker_match) if linker.GetAtomWithIdx(num).GetAtomicNum() != 0 and idx not in f_match]
                if len(set(list(f_match)+list(l_match))) == mol.GetNumHeavyAtoms():
                #if len(set(list(frag_match)+list(linker_match))) == mol.GetNumHeavyAtoms():
                    break
            # Add frag indices
            sub_idx += frag_match
            # Add linker indices to end
            sub_idx += [idx for num, idx in enumerate(linker_match) if linker.GetAtomWithIdx(num).GetAtomicNum() != 0 and idx not in sub_idx]
    
            nodes_to_keep = [i for i in range(len(frag_match))]
    
            aligned_mols.append(Chem.rdmolops.RenumberAtoms(mol, sub_idx))
            aligned_mols.append(frags)
    
            # Renumber dummy atoms to end
            dummy_idx = []
            for atom in aligned_mols[1].GetAtoms():
                if atom.GetAtomicNum() == 0:
                    dummy_idx.append(atom.GetIdx())
            for i, amol in enumerate(aligned_mols):
                sub_idx = list(range(aligned_mols[1].GetNumHeavyAtoms()+2))
                for idx in dummy_idx:
                    sub_idx.remove(idx)
                    sub_idx.append(idx)
                if i == 0:
                    mol_range = list(range(amol.GetNumHeavyAtoms()))
                else:
                    mol_range = list(range(amol.GetNumHeavyAtoms()+2))
                idx_to_add = list(set(mol_range).difference(set(sub_idx)))
                sub_idx.extend(idx_to_add)
                aligned_mols[i] = Chem.rdmolops.RenumberAtoms(amol, sub_idx)
    
            # Get exit vectors
            exit_vectors = []
            linker_atom_idx = []
            for atom in aligned_mols[1].GetAtoms():
                if atom.GetAtomicNum() == 0:
                    if atom.GetIdx() in nodes_to_keep:
                        nodes_to_keep.remove(atom.GetIdx())
                    for nei in atom.GetNeighbors():
                        exit_vectors.append(nei.GetIdx())
                    linker_atom_idx.append(atom.GetIdx())
    
            # Get coords
            conf = aligned_mols[0].GetConformer()
            exit_coords = []
            for exit in exit_vectors:
                exit_coords.append(np.array(conf.GetAtomPosition(exit)))
            linker_coords = []
            for linker_atom in linker_atom_idx:
                linker_coords.append(np.array(conf.GetAtomPosition(linker_atom)))
    
            # Create conformer object for frags
            conf = Chem.Conformer(frags.GetNumAtoms())
            ref_conf = aligned_mols[0].GetConformer()
            for i in range(frags.GetNumAtoms()):
                pos = list(ref_conf.GetAtomPosition(i))
                conf.SetAtomPosition(i, pos)
            conf.SetId(0)
            cid = aligned_mols[1].AddConformer(conf)
    
            # Align frags
            rms = AllChem.AlignMol(aligned_mols[1], aligned_mols[0], atomMap=[(x, x) for x in range(frags.GetNumAtoms())])
            
            edMol = Chem.EditableMol(aligned_mols[1])
            #Remove dummy atoms
            atomsToRemove = []
            for atom in aligned_mols[1].GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atomsToRemove.append(atom.GetIdx())
            
            for a in sorted(atomsToRemove, reverse = True):
                edMol.RemoveAtom(a)
    
            molToUse = edMol.GetMol()
            Chem.SanitizeMol(molToUse)
    
            #return molToUse_Exit, molToUse
            return molToUse
    
        except:
            print(Chem.MolToSmiles(mol), smi_linker, smi_frags)
            return False


    
    ####Hotspots Processing####
    
    
    def processHotspots(self, results, exitVecPos, fragMol, sizeCutOff = 8, minDistance = 1.5, maxDistance = 5, DAthreshold = 10, ApolarThreshold = 1, verbose = False):
        
        #Get density grids
        aboveThresholdDictDonor = results.super_grids['donor'].grid_value_by_coordinates(threshold = DAthreshold)
        aboveThresholdDictAcceptor = results.super_grids['acceptor'].grid_value_by_coordinates(threshold = DAthreshold)
        aboveThresholdDictApolar = results.super_grids['apolar'].grid_value_by_coordinates(threshold = ApolarThreshold)
        
        
        
        if not verbose:
        
            #Check that the fragment lies in a piece of apolar density:
            self.fragInApolar = self.compareFragToApolar(fragMol, aboveThresholdDictApolar)
            
            #Create Donor and Acceptor Mols
            donorMol = self.createDAMols(aboveThresholdDictDonor, aboveThresholdDictApolar, 
                                   exitVecPos, fragMol, sizeCutOff=sizeCutOff, donor = True, 
                                    minDistance=minDistance, maxDistance=maxDistance)
            
            acceptorMol = self.createDAMols(aboveThresholdDictAcceptor, aboveThresholdDictApolar, 
                                   exitVecPos, fragMol, sizeCutOff=sizeCutOff, donor = False, 
                                    minDistance=minDistance, maxDistance=maxDistance)
            
            return {'Donor':donorMol, 'Acceptor':acceptorMol}
    
        elif verbose:
            #Return the intermediate steps of the processing as well
            
            #Check that the fragment lies in a piece of apolar density:
            self.fragInApolar = self.compareFragToApolar(fragMol, aboveThresholdDictApolar)
            
            #Create Donor and Acceptor Mols
            donorOut = self.createDAMols(aboveThresholdDictDonor, aboveThresholdDictApolar, 
                                   exitVecPos, fragMol, sizeCutOff=sizeCutOff, donor = True, 
                                    minDistance=minDistance, maxDistance=maxDistance, verbose = True)
            
            acceptorOut = self.createDAMols(aboveThresholdDictAcceptor, aboveThresholdDictApolar, 
                                   exitVecPos, fragMol, sizeCutOff=sizeCutOff, donor = False, 
                                    minDistance=minDistance, maxDistance=maxDistance, verbose = True)
            
            return donorOut, acceptorOut
        
    
    def createDAMols(self, grid, apolarGrid, exitVecPos, fragMol, sizeCutOff = 8, donor = True, minDistance = 1, maxDistance = 5, verbose = False):
        #Takes a Donor or Acceptor grid and returns an rdkit molecule, where each atom represents a density cluster
        
        #Grid - a grid obtained in processHotspots
        #Donor - True if a Donor grid, False if acceptor
        #apolarGrid - Used to identify if the density clusters are in a favourable position
        #exitVecPosition - 3D coordinates of the exit vector
        #fragMol - a 3D molecule of the fragment 
        #sizeCutOff - The size a density cluster has to be to include it in the molecule output
        #maxDistance - The furthest a donor/acceptor region can be away from the exit vector and still be included
        
        x = []
        y = []
        z = []
    
        for k in list(grid.keys()):
            x.append(grid[k][0][0])
            y.append(grid[k][0][1])
            z.append(grid[k][0][2])
        
        coords = np.transpose(np.array([x,y,z]))
        
        distanceFromExit = []
        for i in coords:
            distanceFromExit.append(self.vectorDistance(exitVecPos, i))
            
        #Only keep density close to the exit vector (but not too close!)
        coordsWithDistance = pd.DataFrame({'x': coords[:,0], 'y' : coords[:,1], 'z' : coords[:,2], 'dist':distanceFromExit})
        within_N_Angstroms = coordsWithDistance.loc[coordsWithDistance['dist'] < maxDistance]
        within_N_Angstroms = within_N_Angstroms.loc[within_N_Angstroms['dist'] > minDistance]
        
        if verbose:
            intermediate1 = within_N_Angstroms #Density after thresholding for distance
        
        #Only keep density such that the fragment atom which is closest to the voxel is the exit point
        within_N_Angstroms = self.compareEVToOtherAtoms(within_N_Angstroms, exitVecPos, fragMol)
        within_N_Angstroms = within_N_Angstroms.loc[within_N_Angstroms['closestToEV']]
        
    
        if verbose:
            intermediate2 = within_N_Angstroms #Density which is closest to ev
    
    
        #Do clustering
        cluster = self.clustering(within_N_Angstroms)
        cluster.clusterDataset()
        cluster.obtainClusterCentroids()
        
        
        centroidsX = [cluster.centroids[key]['center'][0] for key in cluster.centroids.keys()]
        centroidsY = [cluster.centroids[key]['center'][1] for key in cluster.centroids.keys()]
        centroidsZ = [cluster.centroids[key]['center'][2] for key in cluster.centroids.keys()]
        centroidSize = [cluster.centroids[key]['count'] for key in cluster.centroids.keys()]
    
        centroidDF = pd.DataFrame({'x':centroidsX, 'y':centroidsY, 'z':centroidsZ, 'size':centroidSize})
        
        if verbose:
            intermediate3 = centroidDF #Cluster Centroids
        
        
        #Filter hotspots based on their proximity to each other
        #We don't want lots of small hotspots in the same regions as this will massively increase computational time for little benefit
    
        #If centroidDF is empty then just return an empty mol
        
        if centroidDF.shape[0] == 0:
            
            if not verbose:
                return Chem.RWMol()
            else: 
                return Chem.RWMol(), intermediate1, intermediate2, intermediate3
    
    
        elif centroidDF.shape[0] > 0:
            #Otherwise filter further
            centroidDF = self.filterCloseHotspots(centroidDF)
            
            if verbose:
                intermediate4 = centroidDF #Removed cluster centroids which are too close to each other
    
            #If there are no hotspots of the size we want, then check for smaller hotspots
    
            if centroidDF.loc[centroidDF['size'] > sizeCutOff].shape[0] > 0:
                #Only keep centroids which are close to a piece of apolar density (otherwise deemed not to be a promising target)
                centroidDF = self.compareCentroidToApolar(centroidDF, apolarGrid)
                centroidDF = centroidDF.loc[centroidDF['apolarProximity'] < 1] 
            
                if verbose:
                    intermediate5 = centroidDF.loc[centroidDF['size'] > sizeCutOff] #Filtered for apolar proximity
            
            
                if donor:
                    outMol = self.hotspotsDFToMol(centroidDF, cutoff = sizeCutOff, atomType='I') #Iodine atoms represent donor density (arbitrary choice)
                else:
                    outMol = self.hotspotsDFToMol(centroidDF, cutoff = sizeCutOff, atomType='P') #Phosphorus atoms represent acceptor density 
                
                
                if not verbose:
                    return outMol
                else:
                    return outMol, intermediate1, intermediate2, intermediate3, intermediate4, intermediate5
    
            elif centroidDF.loc[centroidDF['size'] > np.floor(sizeCutOff/2)].shape[0] > 0: #Check whether there are smaller hotspot regions if no larger ones
                #Only keep centroids which are close to a piece of apolar density (otherwise deemed not to be a promising target)
                centroidDF = self.compareCentroidToApolar(centroidDF, apolarGrid)
                centroidDF = centroidDF.loc[centroidDF['apolarProximity'] < 1] 
            
                if verbose:
                    intermediate5 = centroidDF.loc[centroidDF['size'] > np.floor(sizeCutOff/2)] #Filtered for apolar proximity
            
                if donor:
                    outMol = self.hotspotsDFToMol(centroidDF, cutoff = np.floor(sizeCutOff/2), atomType='I') #Iodine atoms represent donor density (arbitrary choice)
                else:
                    outMol = self.hotspotsDFToMol(centroidDF, cutoff = np.floor(sizeCutOff/2), atomType='P') #Phosphorus atoms represent acceptor density 
                
                if not verbose:
                    return outMol
                else:
                    return outMol, intermediate1, intermediate2, intermediate3, intermediate4, intermediate5
                
            
            elif centroidDF.loc[centroidDF['size'] > np.floor(sizeCutOff/4)].shape[0] > 0: #Check whether there are smaller hotspot regions if no larger ones
                #Only keep centroids which are close to a piece of apolar density (otherwise deemed not to be a promising target)
                centroidDF = self.compareCentroidToApolar(centroidDF, apolarGrid)
                centroidDF = centroidDF.loc[centroidDF['apolarProximity'] < 1]
    
                if verbose:
                    intermediate5 = centroidDF.loc[centroidDF['size'] > np.floor(sizeCutOff/4)] #Filtered for apolar proximity
    
                if donor:
                    outMol = self.hotspotsDFToMol(centroidDF, cutoff = np.floor(sizeCutOff/4), atomType='I') #Iodine atoms represent donor density (arbitrary choice)
                else:
                    outMol = self.hotspotsDFToMol(centroidDF, cutoff = np.floor(sizeCutOff/4), atomType='P') #Phosphorus atoms represent acceptor density 
    
                if not verbose:
                    return outMol
                else:
                    return outMol, intermediate1, intermediate2, intermediate3, intermediate4, intermediate5
            else:
                #if we're not going to be returning any hotspots
                if not verbose:
                    return Chem.RWMol() # Just return an empty molecule
                else:
                    return Chem.RWMol(), intermediate1, intermediate2, intermediate3, intermediate4
                    



   

    def compareCentroidToApolar(self, centroidDF, apolarGrid):
        
        x = []
        y = []
        z = []
    
        for k in list(apolarGrid.keys()):
            x.append(apolarGrid[k][0][0])
            y.append(apolarGrid[k][0][1])
            z.append(apolarGrid[k][0][2])
        
        coords = np.transpose(np.array([x,y,z]))
        
        proximity = []
        for idx, row in centroidDF.iterrows():
            cPos = np.array([row['x'], row['y'], row['z']])
            
            distanceFromCentroid = []
            for i in coords:
                distanceFromCentroid.append(self.vectorDistance(cPos, i))
            
            proximity.append(min(distanceFromCentroid))
            
        centroidDF['apolarProximity'] = proximity #The nearest bit of apolar density to each centroid
        
        return centroidDF
    
    def compareFragToApolar(self, frag, apolarGrid):
        
        #Return True if all fragment atoms are in close proximity to a piece of apolar density
        #Otherwise return False
        
        x = []
        y = []
        z = []
    
        for k in list(apolarGrid.keys()):
            x.append(apolarGrid[k][0][0])
            y.append(apolarGrid[k][0][1])
            z.append(apolarGrid[k][0][2])
        
        apolar_coords = np.transpose(np.array([x,y,z]))
        frag_coords = np.array([np.array(frag.GetConformer().GetAtomPosition(x)) for x in range(frag.GetNumHeavyAtoms())])
        
        
        #Cross reference apolar_coords against frag_coords:
        

        min_dist_vec = []
        for i in frag_coords:
            
            dist_vec = []
            for j in apolar_coords:
                dist_vec.append(self.vectorDistance(i, j))
            
            min_dist_vec.append(np.min(dist_vec)) #Add the distance between the fragment atom and the closest piece of apolar density
            
        if np.max(min_dist_vec) > 1.5:
            #i.e. there is a fragment atom that is not close to a piece of apolar density
            return False
        else:
            return True
        

    def compareEVToOtherAtoms(self, df, evp, frag):
        
        closestToEV = []
        for idx, row in df.iterrows():
            dPos = np.array([row['x'], row['y'], row['z']])
            distanceToAtoms = []
            for atom in frag.GetAtoms():
                atomPos = np.array(frag.GetConformer().GetAtomPosition(atom.GetIdx()))
                distanceToAtoms.append(self.vectorDistance(dPos, atomPos))
            
            if min(distanceToAtoms) < self.vectorDistance(dPos, evp) - 0.1: #If the distance to an atom is less than the
                                                                        #distance to the exit vector (minus some tolerance)
                closestToEV.append(False)
            else:
                closestToEV.append(True)
            
        df['closestToEV'] = closestToEV
        return df
            

    class clustering:
        
        def __init__(self, df):
            
            self.df = df.reset_index(drop = True)
            
            
            self.df['assignedCluster'] = [-1]*self.df.shape[0]
            
            #Initialise first cluster
            #self.df['assignedCluster', 0] = 0
            self.df.at[0, 'assignedCluster'] = 0

            self.currentClusterIdx = 0
            
        def vectorDistance(self, p1, p2):
            return np.linalg.norm(p1 - p2)
        
        def addToCluster(self, cutoff = 1):
            
            #Get DF of all observations in current cluster
            self.currentClusterDF = self.df.loc[self.df['assignedCluster'] == self.currentClusterIdx][['x', 'y', 'z']]
            
            distVec = []
            for idx, row in self.currentClusterDF.iterrows():
                
                #Compute distance to all unassigned datapoints
                for idx2, row2 in self.df.loc[self.df['assignedCluster'] == -1][['x', 'y', 'z']].iterrows():
                    dist = self.vectorDistance(np.array([row2['x'], row2['y'], row2['z']]), np.array([row['x'], row['y'], row['z']]))
                    distVec.append(dist)
                    
                    if dist <= cutoff:
                        #We've found an observation that's close to the current cluster
                        #self.df['assignedCluster'][idx2] = self.currentClusterIdx
                        self.df.at[idx2, 'assignedCluster'] = self.currentClusterIdx
                        return 0
            
            #If we haven't found an observation that's close then we update the current cluster by one and assign a first value
            if min(distVec) > cutoff:
                self.currentClusterIdx += 1
                
                for idx, row in self.df.loc[self.df['assignedCluster'] == -1].iterrows():
                    #self.df['assignedCluster'][idx] = self.currentClusterIdx
                    self.df.at[idx, 'assignedCluster'] = self.currentClusterIdx
                    break
        
        
        def clusterDataset(self):
            
            for _ in range(self.df.shape[0] - 1):
                self.addToCluster()
                
        def obtainClusterCentroids(self):
            
            self.centroids = {}
            
            for idx in self.df['assignedCluster'].drop_duplicates():
                self.centroids[idx] = {'center':list(self.df.loc[self.df['assignedCluster'] == idx][['x', 'y', 'z']].mean()), 
                                       'count':self.df.loc[self.df['assignedCluster'] == idx].shape[0]}
            
        

    def hotspotsDFToMol(self, df, cutoff = 0, atomType = 'I'):
        
        #Read in a cluster centroids df with coordinates for the center and the size
        #Return a molecule with a conformer with all the centroids above a certain size
            #(default size = 0)
            
        df = df.loc[df['size'] > cutoff].reset_index(drop = True)
        #Create the molecule with the correct number of atoms
        fakelig = atomType
        for i in range(df.shape[0] - 1):
            fakelig += f'.{atomType}'
        molPharmProf = Chem.MolFromSmiles(fakelig)
        
        #Now add the 3d position to each atom
        conf = Chem.Conformer(molPharmProf.GetNumAtoms())
        
        for idx, row in df.iterrows():
            conf.SetAtomPosition(idx, list(row[['x', 'y', 'z']]))
    
        conf.SetId(0)
        molPharmProf.AddConformer(conf)
        
        return molPharmProf

    def vectorDistance(self, p1, p2):
        return np.linalg.norm(p1 - p2)

    def vectorAngle(self, x, y):
        #Returns the angle between two numpy arrays
        cos_theta = np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))
    
        return np.arccos(cos_theta)

    def filterCloseHotspots(self, df, distThreshold = 1.5):
        #Do filtering of hotspots based on their proximity to each other
        #Take the largest hotspot and add to list
        #Compute distance from that hotspot to all other hotspots and add to the list the largest hotspot with a greater distance than the distance threshold
        #Compute distance of all remaining hotspots to hotspots in the list and add to the list the largest hotspot with a greater distance than the distance threshold to all hotspots currently in the list.
        
        df = df.sort_values('size', ascending = False) #sort based on size of hotspot
        
        filtered = pd.DataFrame()
        filtered = filtered.append(df.loc[df['size'] == max(df['size'])]) #Add largest hotspot
        
        for idx, row in df.iterrows():
        
            add = 1
            for idx2, row2 in filtered.iterrows():
                if self.vectorDistance(np.array([row['x'], row['y'], row['z']]), np.array([row2['x'], row2['y'], row2['z']])) < distThreshold:
                    add = 0 #The hotspot under consideration is too close to one that's already included in the list
    
            if add == 1: #i.e. the hotspot isn't too close to any of the ones already in the filtered list, then we add it to the filtered list
                filtered = filtered.append(row)
                
                
        return filtered
    
    
    def prepareProfiles(self, HotspotsDF, single = True):

        #Set up the counts profiles corresponding to the hotspots
        #We'll use those profiles to make elaborations

        if single:

            hSingles = {}
            for i, r in HotspotsDF.iterrows():
                hSingles[f'{i}'] = {'type':r['type'], 'position':r['position'], 'distFromExit':r['distFromExit'], 'angFromExit':r['angFromExit']}
                
            return hSingles
        
        else:
            #Prepare profile for satisfying multiple pharmacophores
        
            HotspotsDF = HotspotsDF.sort_values('distFromExit', ascending = False, ignore_index = True) 
            hMulti = {}
            
            for i, r in HotspotsDF.iterrows():
                hMulti[i] = {'type':r['type'], 'position':r['position'], 'distFromExit':r['distFromExit'], 'angFromExit':r['angFromExit']}
        
            return hMulti
        
        '''
        self.hPairs = {}
        for i, r in self.HotspotsDF.iterrows():
            for ii, rr in self.HotspotsDF.iterrows():

                if r['distFromExit'] > rr['distFromExit']:
                    self.hPairs[f'{ii}To{i}'] ={1:{'type': rr['type'], 'position': rr['position'], 'distFromExit':rr['distFromExit'], 'angFromExit':rr['angFromExit']}, 2:{'type': r['type'], 'position':r['position'], 'distFromExit':r['distFromExit'], 'angFromExit':r['angFromExit']}}

        '''

        
    
    def getLattice(self, evp, size = 5, resolution = 1):
        x = np.arange(evp[0], evp[0] + size, resolution) - size/2
        y = np.arange(evp[1], evp[1] + size, resolution) - size/2
        z = np.arange(evp[2], evp[2] + size, resolution) - size/2
    
        lattice = np.vstack(np.meshgrid(x,y,z)).reshape(3,-1).T
    
        return lattice




    
    
    
