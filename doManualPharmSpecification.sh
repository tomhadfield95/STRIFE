
mol3DPath=$1
fragSmilesLoc=$2
proteinLoc=$3
outStoreDir=$4

echo 'Loading PyMol for Manual Pharmacophore Specification'



mkdir ${outStoreDir}

source <path/to/CSD/Suite>/Python_API_2021/miniconda/bin/activate STRIFE
python prepareLattice.py ${mol3DPath} ${fragSmilesLoc} ${proteinLoc} ${outStoreDir}



source <path/to/CSD/Suite>/Python_API_2021/miniconda/bin/activate pymol_env

python loadLatticeIntoPymol.py ${outStoreDir}




