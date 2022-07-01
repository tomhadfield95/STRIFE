#Usage:

#IMPORTANT: You must change the penultimate line - replace <path/to/conda/installation> with the location of your miniconda or anaconda installation
	#   You can run conda info --envs to find where your conda environments are installed

#Alternatively you can activate a different conda environment which has a working version of PyMol installed


#run 'bash do_manual_pharm_specification <path/to/Mol/SDF> <smiles string> <path/to/protein/pdb> <output/directory>'

#where: <path/to/mol/SDF> is an SDF of either the fragment you want to elaborate from, or a larger molecule of which the desired fragment is a substructure
       #<smiles string> can either be a file containing the SMILES of the desired fragment or a string containing the fragment - the fragment exit vector should be denoted by a dummy atom
       #<path/to/protein/pdb> is the path to the protein pdb file
       #<output/directory> is a directory which files needed to run the PyMol session will be stored

#Once the PyMol session has started - select the points in the lattice that you wish to use as pharmacophoric points and save them in the output directory. Donor points must be saved as an SDF molecule called
	#donorHotspot.sdf and acceptor points must be saved as acceptorHotspot.sdf

#WHEN SAVING THE SELECTED IN LATTICE POINTS IN PyMol, YOU MUST USE THE state=0 command - i.e. (in the PyMol command line): save donorHotspot.sdf, my_HBD_selection, state=0

#Once the SDF files have been saved, you can use them to run STRIFE - to do so you must specify --model_type 1 and --output_directory <output/directory> (i.e. the same place as you saved everything above) as arguments to STRIFE.py




path_to_mol_SDF=$1
fragment_smiles=$2
path_to_protein_pdb=$3
output_directory=$4

echo 'Loading PyMol for Manual Pharmacophore Specification'



mkdir ${output_directory}

#source <path/to/conda/installation>/bin/activate STRIFE #MUST CHANGE THIS FILE PATH TO MATCH YOUR OWN INSTALLATION!!!
python prepareLattice.py ${path_to_mol_SDF} ${fragment_smiles} ${path_to_protein_pdb} ${output_directory}


#source <path/to/conda/installation>/bin/activate pymol_env #MUST CHANGE THIS FILE PATH TO MATCH YOUR OWN INSTALLATION!!!
python loadLatticeIntoPymol.py ${output_directory}




