#A file to store functions for making elaborations based on a specific profile

#Import Libaries

#########Standard Libraries##########
import numpy as np
import pandas as pd



#########RDKit Modules############
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMMPA
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') #Supress annoying RDKit output
from multiprocessing import Pool

##########
import frag_utils





savedModelDict = {'Orig':'GenModelDeLinker_saved.pickle',
                 'Count':'coarseGrainedGenModel_saved.pickle',
                 'Pharm':'fineGrainedGenModel_saved.pickle'}



def valueCountsToMols(valueCounts, howmany = None):

    if howmany is None:
        howmany = valueCounts.shape[0] #if we don't specify a subset size, we'll take them all

    mols = list(pd.Series(list(valueCounts.index)).apply(Chem.MolFromSmiles)) #Pulls the index, which is the molecule associated with the value count
    img = Draw.MolsToGridImage(mols[:howmany], legends = [str(count) for count in valueCounts]) #Creates grid with the molecule image and the number of times it was generated

    return img

def remove_dummy_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol2 = AllChem.ReplaceSubstructs(mol, Chem.MolFromSmiles('*'), Chem.MolFromSmiles('[H]'), True)[0]
    mol3 = Chem.RemoveHs(mol2)

    return Chem.MolToSmiles(mol3)

def conEmbed(gen, fragConfFile):

    refConf = Chem.SDMolSupplier(fragConfFile)[0]

    genMols = []
    for i in range(50):
        genMol = Chem.MolFromSmiles(gen)
        AllChem.ConstrainedEmbed(genMol, refConf, randomseed=np.random.randint(0, 10000))
        genMols.append(genMol)


    #Find the one with lowest energy:
    energies = [AllChem.UFFGetMoleculeForceField(m).CalcEnergy() for m in genMols]

    molToUse = genMols[np.argmin(energies)]

    return molToUse



def getValueCountsImage(data, substructMatchMol = None):

    if substructMatchMol is not None:
        data['subs'] = [Chem.MolFromSmiles(s).HasSubstructMatch(substructMatchMol) for s in data['gen']]

        return valueCountsToMols(data.loc[data['subs']]['gen'].value_counts())

    else:

        return valueCountsToMols(data['gen'].value_counts())





def getMultipleSDFs(data, core, outNameCore):

    #Iterate over the different indices (i.e. the different pharmacophoric profiles) and pass to getSDFs
    
    for idx in data['idx'].drop_duplicates():
        
        d = data.loc[data['idx'] == idx] 
        
        getSDFs(d, core, f'{outNameCore}_p{idx}.sdf')
    

def getSDFs(data, core, outName):
    #Core is the embedded mol we use for the constrained embedding
    
    mols = []
    dfIdx = []
    sdfIdx = []
    sdfIdxNum = 0
    for idx, g in enumerate(list(data['gen'])):
        dfIdx.append(idx)
        m = Chem.MolFromSmiles(g)
#         AllChem.ConstrainedEmbed(m, core)
        
        try:
            AllChem.ConstrainedEmbed(m, core)
            fail = 0
        except:
            fail = 1
            print(idx)


        #Now addHs and reembed:
        mHs = Chem.AddHs(m)
#         AllChem.ConstrainedEmbed(mHs, m)
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


#Do Filtering of generated molecules

def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def addElab(df):

    elab = []
    for idx, row in df.iterrows():
        try:
            elab.append(returnGenElaborationWithDummy(row['frag'], row['gen']))
        except:
            elab.append('Error')


    df['elab'] = elab

    return df




def filterGeneratedMols(df):

    df = parallelize_dataframe(df, addElab, n_cores = 7)
    df = frag_utils.check_2d_filters_dataset_pandas(df)

    df = df.loc[df['passed'] == 1][['frag', 'full', 'gen', 'elab']]

    return df


def vector_angle(x,y):
    #Returns the angle between two numpy arrays
    cos_theta = np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))

    return np.arccos(cos_theta)

def vector_distance(x,y):
    #Returns the distance between two numpy arrays
    diff = np.subtract(x,y)
    return np.linalg.norm(diff)






def returnGenElaboration(frag, gen):
    #Both inputs are smiles strings - frag should only have a single dummy atom
    #Return the elaboration
    #i.e. the part of gen which is not part of frag

    generated_elaboration = single_elaboration_2(frag, gen)

    if generated_elaboration == 'Unable to determine elaboration':
        generated_elaboration = single_elaboration(frag, gen)

    #Both molecules have dummy atoms which we want to remove:


    try:
        generated_without_dummy = remove_dummy_atom(generated_elaboration)
    except:
        with open('elaborations_debug.txt', 'w+') as f:
            f.write(generated_elaboration)
            generated_without_dummy = 'Unable to compute - placeholder string'

    return generated_without_dummy

def returnGenElaborationWithDummy(frag, gen):
    #Both inputs are smiles strings - frag should only have a single dummy atom
    #Return the elaboration
    #i.e. the part of gen which is not part of frag

    generated_elaboration = single_elaboration_2(frag, gen)

    if generated_elaboration == 'Unable to determine elaboration':
        generated_elaboration = single_elaboration(frag, gen)

    return generated_elaboration


def single_elaboration_2(frag, full):
    #take a molecule and fragment of that molecule - return the elaborated structure by enumerating cuts
    #to the molecule and working with the pair that matches the fragment

    output_smiles = None

    cuts = fragment_mol_for_elaboration(full, full)

    #Now iterate through:

    for _, _, frag1, frag2 in cuts:
        frag1, frag2
        c1 = compare_by_inchi(frag1, frag)
        c2 = compare_by_inchi(frag2, frag)

        if c1 == 1:
            output_smiles = frag2
        if c2 == 1:
            output_smiles = frag1

    if output_smiles is None:
        return 'Unable to determine elaboration'


    return output_smiles

def single_elaboration(frag, full):
    #Aligns by fragments, finds the exit node and fragments on that.
    #Returns the molecule which doesn't have a substruct match with frag


    #FRAG SHOULD ONLY HAVE A SINGLE DUMMY ATOM
    if len(frag.split('.')) > 1:
        frag = frag.split('.')[0]


    length_of_elab = Chem.MolFromSmiles(full).GetNumHeavyAtoms() - Chem.MolFromSmiles(frag).GetNumHeavyAtoms()

    if length_of_elab < 2:
        return '*'
    else:

        (aligned_mol, aligned_frag), to_keep, exit = align_smiles_by_frags(full, frag)

        #Get neighbours of exit atom
        neighbors = aligned_mol.GetAtomWithIdx(exit[0]).GetNeighbors()

        #Iterate through neighbours and take the one which is in the elaborated section
        for neighbor in neighbors:
            if neighbor.GetIdx() not in to_keep:
                bond = aligned_mol.GetBondBetweenAtoms(aligned_mol.GetAtomWithIdx(exit[0]).GetIdx(), neighbor.GetIdx())

        new_mol = Chem.FragmentOnBonds(aligned_mol, [bond.GetIdx()])

        #We want to make sure we take the elaborated bit
        # Include dummy in query
        du = Chem.MolFromSmiles('*')
        qp = Chem.AdjustQueryParameters()
        qp.makeDummiesQueries=True
        qfrag = Chem.AdjustQueryProperties(Chem.MolFromSmiles(frag),qp)



        #new_mol.HasSubstructMatch(qfrag)
        new_mol_smiles = Chem.MolToSmiles(new_mol).split('.')
        #return new_mol_smiles

        #Check the two fragments
        finished = 0
        for smiles in new_mol_smiles:
            if Chem.MolFromSmiles(smiles).HasSubstructMatch(qfrag) == False:
                return smiles
            else:
                print(smiles)


        for smiles in new_mol_smiles:
            #Get rid of dummy atoms
            smiles2 = remove_dummy_atom(smiles)
            #Compute inchikey
            inchisplit = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(smiles2))


            #Now get the inchi key of the original fragment
            frag2 = remove_dummy_atom(frag)
            inchifrag = Chem.inchi.MolToInchiKey(Chem.MolFromSmiles(frag2))

            if inchisplit != inchifrag:
                return smiles

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

def fragment_mol_for_elaboration(smi, cid, pattern="[#6+0;!$(*=,#[!#6])]!@!=!#[*]"):
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

def compare_by_inchi(mol1, mol2):
    #take two smiles strings and assess if they equal by converting to inchi strings
    #helper function for single_elaboration_2
    #Remove dummy atoms
    mol1 = mol1.split('.')[0]
    mol2 = mol2.split('.')[0]

    mol1 = remove_dummy_atom(mol1)
    mol2 = remove_dummy_atom(mol2)



    inchi1 = Chem.MolToInchi(Chem.MolFromSmiles(mol1))
    inchi2 = Chem.MolToInchi(Chem.MolFromSmiles(mol2))


    if inchi1 == inchi2:
        return 1
    else:
        return 0

def parallelize_dataframe(df, func, n_cores=4):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def addElab(df):

    elab = []
    for idx, row in df.iterrows():
        try:
            elab.append(returnGenElaborationWithDummy(row['frag'], row['gen']))
        except:
            elab.append('Error')


    df['elab'] = elab

    return df


def checkIfNovel(df, db):
    
    novel = []
    
    for idx, row in df.iterrows():

        if row['elab'] == 'Error' or row['elab'] == 'Unable to compute - placeholder string':
            novel.append(-1)
        elif remove_dummy_atom(row['elab']) in db:
            novel.append(0)
        else:
            novel.append(1)
           
    df['novel'] = novel
    return df


def checkIfValid(df):

    valid = []

    for idx, row in df.iterrows():
        if Chem.MolFromSmiles(row['frag']).GetNumHeavyAtoms() < Chem.MolFromSmiles(row['gen']).GetNumHeavyAtoms():
            valid.append(1)
        else:
            valid.append(0)
    
    df['valid'] = valid

    return df



'''
def checkProfileCount(mol, frag, profile):
    #Take a molecule as input and check whether it has the same number of pharmacophores as specified in the profile

    #Profile should be of the form: [dist (=0), ang (=0), HBANum, HBDNum, AromaticNum]
    
    molPharms = getPharmNumbers(mol)
    fragPharms = getPharmNumbers(frag)
    
    elabProfile = [0, 0, molPharms['Acceptor'] - fragPharms['Acceptor'], molPharms['Donor'] - fragPharms['Donor'], molPharms['Aromatic'] - fragPharms['Aromatic']]
    
    if elabProfile == profile:
        return True
    else:
        return False
    
    
    
def getPharmNumbers(mol):
    
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

def replacePharmProfile(newProfile, jsonF, replacementType='Orig'):

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


def generateMolsFromSpecificProfile(frag, elabLength, modelType, profile, smiName, savedModelDict = savedModelDict):

    #First want to create a model outside of the function so that we only need to update the model.valid_data attribute
    #Instead of loading the entire model repeatedly

    smi = createSmi(frag, elabLength)
    smi.to_csv(f'{smiName}.smi', header = False, index = False, sep = ' ')

    smi.columns = ['full', 'chopped', 'frag', 'dist', 'ang']
    smi['pharm'] = [addPharmProfile.createPharmacophoricProfile(row['frag'], row['full']) for idx, row in smi.iterrows()]
    smi.to_csv(f'{smiName}Pharm.smi', header = None, index = None, sep = ' ')

    raw_data = train_valid_split(f'{smiName}Pharm.smi')['valid']

    jsonFile = preprocessSingle(raw_data)

    #Replace pharmacophoric profile
    jsonFile = replacePharmProfile(profile, jsonFile, modelType)



    jsonName = f'{smiName}Json{modelType}.json'


    with open(jsonName, 'w') as f:
        json.dump(jsonFile, f)

    data_input = jsonName



    with open(data_input, 'r') as f:
        data = json.load(f)[0]

    config_input = '{"generation":true, "number_of_generation_per_valid": 250, "num_timesteps": 9, "hidden_size": 50, "encoding_size": 4, "batch_size": 1, "qed_trade_off_lambda": 0, "node_keep_trade_off_lambda": 0, "fwd_bkd_trade_off_lambda": 0, "num_epochs": 2, "epoch_to_generate":2, "compensate_num": 0, "num_different_starting": 1, "train_file":' + '"' + data_input + '"' + ', "valid_file":' + '"' + data_input + '"' + ', "output_dir": "store_out_temp", "random_seed": 4}'
    args = {'--dataset':'zinc', '--config':config_input,  '--restore': savedModelDict[modelType]}
    dataset = args.get('--dataset')



    if modelType == 'Orig':
        model = DenseGGNNChemModelOrig(args)
    elif modelType == 'Count':
        model = DenseGGNNChemModelCount(args)
    elif modelType == 'Pharm':
        model = DenseGGNNChemModelPharm(args)


    try:
        os.remove("generatedMols.smi") #Delete molecules which generated in a previous run
    except:
        pass
    #Generate Molecules

    model.train()

    outDF = pd.read_csv(f"{model.run_id}_generated_smiles_{model.params['dataset']}.smi", header = None, sep = ' ')
    outDF.columns = ['frag', 'full', 'gen']
    return outDF


def addAlkane(fragment, alkLength):

    #Create alkane to add to the fragment
    alk = 'C'*int(alkLength) + '[*:1]'
#     alk = 'C'*3 + '[*:1]'

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


def createSmi(fragment, alkLength):

    full = addAlkane(fragment, alkLength)
    chopped = 'C'*int(alkLength) + '[*:1]'
#     chopped = 'C'*3 + '[*:1]'
    smi = pd.DataFrame({'full':full, 'chopped':chopped, 'frag':fragment, 'dist':[0], 'ang':[0]})

    return smi

def preprocessSingle(raw_data):
    processed_data = []

    file_count = 0

    for i, (smiles_mol, smiles_frag, QED_in, abs_dist) in enumerate([(mol['smiles_mol'], mol['smiles_frag'],
                                                                           mol['QED_in'], mol['abs_dist']) for mol in raw_data]):
        #print(file_count, end='\r')
        #print("File %d" % i)
        (mol_out, mol_in), nodes_to_keep, exit_points = align_smiles_by_frags(smiles_mol, smiles_frag)
        if mol_out == []:
            continue
        nodes_in, edges_in = to_graph_mol(mol_in, dataset)
        nodes_out, edges_out = to_graph_mol(mol_out, dataset)
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



'''




