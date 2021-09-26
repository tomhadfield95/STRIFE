#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:12:51 2021

@author: hadfield
"""

#Implementation of the preprocessing class without any of the features 
#dependant on the CSD package


# Generic imports
import numpy as np
import sys
from itertools import product


#rdkit imports
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMMPA


class preprocessing:
    
    #Class for various preprocessing tasks associated with the STRIFE algorithm.
    
    def __init__(self):
        pass
    
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


    
    
    def vectorDistance(self, p1, p2):
        return np.linalg.norm(p1 - p2)

    def vectorAngle(self, x, y):
        #Returns the angle between two numpy arrays
        cos_theta = np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))
    
        return np.arccos(cos_theta)

    
    
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
        
       
        
    
    def getLattice(self, evp, size = 5, resolution = 1):
        x = np.arange(evp[0], evp[0] + size, resolution) - size/2
        y = np.arange(evp[1], evp[1] + size, resolution) - size/2
        z = np.arange(evp[2], evp[2] + size, resolution) - size/2
    
        lattice = np.vstack(np.meshgrid(x,y,z)).reshape(3,-1).T
    
        return lattice




    
    
    
