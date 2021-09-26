from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def align_smiles_by_MCS_repeat(smiles_1, smiles_2):
    mols = [Chem.MolFromSmiles(smiles_1), Chem.MolFromSmiles(smiles_2)]
    
    # Align mols by MCS
    res=rdFMCS.FindMCS(mols, ringMatchesRingOnly=True)
    
    aligned_mols = []
    partial_mols = []
    for mol in mols:
        sub_idx = list(mol.GetSubstructMatch(Chem.MolFromSmarts(res.smartsString)))
        nodes_to_keep = [i for i in range(len(sub_idx))]
        size_MCS = len(sub_idx)
        #print(size_MCS)
        mol_range = list(range(mol.GetNumHeavyAtoms()))
        idx_to_add = list(set(mol_range).difference(set(sub_idx)))
        sub_idx.extend(idx_to_add)
        aligned_mols.append(Chem.rdmolops.RenumberAtoms(mol, sub_idx))
        partial_mols.append(Chem.rdmolops.DeleteSubstructs(aligned_mols[-1], Chem.MolFromSmarts(res.smartsString)))
        #partial_mols.append(Chem.ReplaceCore(aligned_mols[-1], Chem.MolFromSmarts(res.smartsString), aligned_mols[-1].GetSubstructMatch(Chem.MolFromSmarts(res.smartsString)), labelByIndex=True))
        
    # Align any remaining rings
    # Check if any rings not within MCS
    aligned_mols_3 = []
    ssr_refs = []
    for mol in aligned_mols:
        ssr = Chem.GetSymmSSSR(mol)
        ssr_ref = []
        used_idx = list(range(size_MCS))
        for r in ssr:
            if not bool(set(list(r)).intersection(set(used_idx))):
                ssr_ref.append(list(r))
        ssr_refs.append(ssr_ref)
    # If one ring each not aligned, align it
    size_ring = 0
    if len(ssr_refs[0]) == 1 and len(ssr_refs[1]) == 1:
        # Don't align 3 member ring with 6 member ring
        if abs(len(ssr_refs[0][0]) - len(ssr_refs[1][0])) < 3:
            size_ring = min(len(ssr_refs[0][0]), len(ssr_refs[1][0]))
            for idx, mol in enumerate(aligned_mols):
                sub_idx = list(ssr_refs[idx][0])
                #size_ring = len(sub_idx)
                #print(size_ring)
                mol_range = list(range(mol.GetNumHeavyAtoms()))
                sub_idx = list(range(size_MCS)) + sub_idx
                idx_to_add = list(set(mol_range).difference(set(sub_idx)))
                sub_idx.extend(idx_to_add)
                aligned_mols_3.append(Chem.rdmolops.RenumberAtoms(mol, sub_idx)) 
        else:
            aligned_mols_3 = aligned_mols
    else:
        aligned_mols_3 = aligned_mols
        
    # Find MCS of non-aligned portions of mols
    res_2 = rdFMCS.FindMCS(partial_mols, ringMatchesRingOnly=True)    
    #print("MCS_2 smarts: ", res_2.smartsString)
    # Align mols by res_2
    aligned_mols_2 = []
    for mol in aligned_mols_3:
        # Find match(es) to MCS
        sub_idx = list(mol.GetSubstructMatches(Chem.MolFromSmarts(res_2.smartsString)))
        #print("MCS_2 sub_idx: ", sub_idx)
        # Check not C-C only - if so skip
        #print(Chem.MolToSmiles(Chem.MolFromSmarts(res_2.smartsString)))
        if Chem.MolToSmiles(Chem.MolFromSmarts(res_2.smartsString)) in ['CC', 'cc', 'C:C']:
            sub_idx = []
        # Check if any matches - if not skip
        if sub_idx == [()]:
            sub_idx = []
        # Remove atoms in aligned ring from match
        for idx in range(len(sub_idx)):
            try: 
                #print("Old: ", sub_idx[idx])
                #if size_ring != 0:
                new_sub_idx = []
                for atom_idx in sub_idx[idx]:
                    if atom_idx not in list(range(size_MCS,size_MCS+size_ring)):
                        new_sub_idx.append(atom_idx) 
                    #sub_idx[idx] = set(list(sub_idx[idx])).difference(set(list(range(size_MCS,size_MCS+size_ring))))
                sub_idx[idx] = new_sub_idx
                #print("New: ", sub_idx[idx])
                min(sub_idx[idx]) # If fail then empty set so set =[0]
            except:
                sub_idx[idx] = [0]
        # If multiple matches, check if any overlap with previously aligned mol
        if len(sub_idx) > 1:
            for match in sub_idx:
                #print("Match: ", match)
                if min(match) < (size_MCS+size_ring):
                    continue
                else:
                    sub_idx = list(match)
        elif len(sub_idx) == 1 and min(sub_idx[0]) >= (size_MCS+size_ring):
            #print(sub_idx[0])
            sub_idx = list(sub_idx[0])
        if sub_idx != []:
            if type(sub_idx[0]) != int:
                sub_idx = []
        size_MCS_2 = len(sub_idx)
        mol_range = list(range(mol.GetNumHeavyAtoms()))
        sub_idx = list(range(size_MCS+size_ring)) + sub_idx
        idx_to_add = list(set(mol_range).difference(set(sub_idx)))
        sub_idx.extend(idx_to_add)
        aligned_mols_2.append(Chem.rdmolops.RenumberAtoms(mol, sub_idx))
        
    # Loop over remaining atoms
    # First find atoms coming off the same sidechain
    size_side_chains = 0
    aligned_mols_4 = aligned_mols_2
    #print(size_MCS, size_ring, size_MCS_2, size_side_chains)
    while True:
        side_chain_dict = [{}, {}]
        # Find atoms attached to core
        for idx, mol in enumerate(aligned_mols_4):
            for atom in mol.GetAtoms():
                if atom.GetIdx() < size_MCS + size_ring + size_MCS_2 + size_side_chains:
                    continue
                else:
                    for bond in atom.GetBonds():
                        if bond.GetBeginAtomIdx() != atom.GetIdx() and bond.GetBeginAtomIdx() < size_MCS + size_ring + size_MCS_2 + size_side_chains:
                            if atom.GetIdx() in side_chain_dict[idx]:
                                side_chain_dict[idx][atom.GetIdx()] = min(side_chain_dict[idx][atom.GetIdx()], bond.GetBeginAtomIdx())
                            else:
                                side_chain_dict[idx][atom.GetIdx()] = bond.GetBeginAtomIdx()
                        elif bond.GetEndAtomIdx() != atom.GetIdx() and bond.GetEndAtomIdx() < size_MCS + size_ring + size_MCS_2 + size_side_chains:
                            if atom.GetIdx() in side_chain_dict[idx]:
                                side_chain_dict[idx][atom.GetIdx()] = min(side_chain_dict[idx][atom.GetIdx()], bond.GetEndAtomIdx())
                            else:
                                side_chain_dict[idx][atom.GetIdx()] = bond.GetEndAtomIdx()
        # Check for overlap or break if no overlap
        #print(side_chain_dict[0].values(), side_chain_dict[1].values())
        intersection = set(side_chain_dict[0].values()).intersection(set(side_chain_dict[1].values()))
        #print(list(side_chain_dict[0].values()), list(side_chain_dict[1].values()))
        #print("Side chain intersection: ", intersection)
        if not intersection:
            break
        else:
            sub_idxs = [[], []]
            for val in intersection:
                sub_idxs[0].append(list(side_chain_dict[0].keys())[list(side_chain_dict[0].values()).index(val)])
                sub_idxs[1].append(list(side_chain_dict[1].keys())[list(side_chain_dict[1].values()).index(val)])
            #print(sub_idxs)
        # Align mols and repeat
        for i, mol in enumerate(aligned_mols_2):
            mol_range = list(range(mol.GetNumHeavyAtoms()))
            sub_idx = list(range(size_MCS+size_ring+size_MCS_2+size_side_chains)) + sub_idxs[i]
            #print(sub_idx)
            idx_to_add = list(set(mol_range).difference(set(sub_idx)))
            sub_idx.extend(idx_to_add)
            #print(sub_idx)
            aligned_mols_4[i] = Chem.rdmolops.RenumberAtoms(mol, sub_idx)
        size_side_chains += len(intersection)
        #print("Aligned side chains: ", size_side_chains)
        #break
        
    # Second match remaining atoms by type
    size_matched_types = 0
    while True:
        atom_type_dict = [{}, {}]
        # Find remaining atoms types
        for idx, mol in enumerate(aligned_mols_4):
            for atom in mol.GetAtoms():
                # Check if atom has been aligned
                if atom.GetIdx() < size_MCS + size_ring + size_MCS_2 + size_side_chains + size_matched_types:
                    continue
                # If not, add atom type to dictionary 
                else:
                    atomic_num = atom.GetAtomicNum()
                    # Don't align carbons
                    if atomic_num == 6:
                        continue
                    if atomic_num not in atom_type_dict[idx]:
                        atom_type_dict[idx][atomic_num] = [atom.GetIdx()]
                    else:
                        atom_type_dict[idx][atomic_num].append(atom.GetIdx())
        # Check for atom type overlap or break if no overlap
        intersection = set(atom_type_dict[0].keys()).intersection(set(atom_type_dict[1].keys()))
        #print(atom_type_dict[0].keys(), atom_type_dict[1].keys())
        #print("Atom type intersection: ", intersection)
        if not intersection:
            break
        else:
            sub_idxs = [[], []]
            total_num = 0
            for key in intersection:
                #print(atom_type_dict[0][key], atom_type_dict[1][key])
                num = min(len(atom_type_dict[0][key]), len(atom_type_dict[1][key]))
                total_num += num
                sub_idxs[0].extend(atom_type_dict[0][key][0:num])
                sub_idxs[1].extend(atom_type_dict[1][key][0:num])
            #print(sub_idxs)
        # Align mols and repeat
        for i, mol in enumerate(aligned_mols_2):
            mol_range = list(range(mol.GetNumHeavyAtoms()))
            sub_idx = list(range(size_MCS+size_ring+size_MCS_2+size_side_chains+size_matched_types)) + sub_idxs[i]
            #print(sub_idx)
            idx_to_add = list(set(mol_range).difference(set(sub_idx)))
            sub_idx.extend(idx_to_add)
            #print(sub_idx)
            aligned_mols_4[i] = Chem.rdmolops.RenumberAtoms(mol, sub_idx)
        size_matched_types += total_num
        #print("Aligned atom types: ", size_matched_types)
    
    return aligned_mols_4, res, nodes_to_keep
