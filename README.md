# STRIFE (STRucture Informed Fragment Elaboration)

This repository contains our implementation of the STRIFE algorithm, described in this biorxiv preprint (TODO: Provide link when up).



# Acknowledgements

The generative models in this work were based on the [Constrained Graph Variational Autoencoder](https://arxiv.org/abs/1805.09076). We thank the authors for releasing the [code](https://github.com/Microsoft/constrained-graph-variational-autoencoder) and invite you to cite their work if you found it useful. 

STRIFE calculates Fragment Hotspot Maps using the [Hotspots API](https://github.com/prcurran/hotspots) developed by [Curran et al.](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00996). We thank the authors for making this resource available.



# Requirements

This code was tested in Python 3.7, using Tensorflow 1.15. A yaml file containing all required libraries is provided. Please see 'Installation' below for more information.

Use of the Hotspots API (extraction of target-specific information) and the GOLD API (constrained docking) requires the commercial CSD Python API,  maintained by the [Cambridge Crystallographic Data Centre](https://www.ccdc.cam.ac.uk/) (CCDC). If you do not have a CSD licence, STRIFE allows you to easily specify your own pharmacophoric points (see below) but you will need to provide your own function for docking and scoring the generated molecules.


# Installation

## Downloading the CSD Python API

To download the CSD suite you need to know your institution's Activation Key and Licence Number. Go to the [CSD Downloads page](https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/) and fill in the information to receive an email (titled 'CCDC Download Request') with a set of download links. Copy the relevant link (for linux users 'CSDS 2021.2 Linux' or similar).

As the link is very long, it may be necessary to enclose the wget command in a .sh script, with the link address enclosed by single quotes (see below)

```
wget '<insert download link here>'
```

Then run the .sh script from the command line (assuming you named the bash script `download_CSD.sh`):

```
bash download_CSD.sh
```

If successful, your current directory will now contain a tar file, which can be unzipped using the `tar -zxvf <tar_file>` functionality.

This should yield a directory named `csds-2021.2.0-linux-installer` (or similar, depending on the exact version). You can run the installer by entering `csds-2021.2.0-linux-installer/csds-linux-x64.run` in the command line. If you don't have access to a GUI, you should add the `--mode text` flag.

Hopefully, after this stage you should have access to the CSD suite!



## Downloading GHECOM

The Hotspots API uses GHECOM for buriedness calculations. It can be downloaded [here](https://pdbj.org/ghecom/download_src.html). After providing some basic contact information, you will be provided with a link to download the source code (at time of writing labelled `ghecom-src-20200721.tar.gz`). To install GHECOM, download the link and run the following commands:

```
tar -zxvf ghecom-src-[date].tar.gz
cd src
make
```

If successful, an executable file called `ghecom` will now be in the parent directory.


## Setting up the STRIFE conda environment

This can be done by simply using the provided `STRIFE_environment.yml` file. You will need to edit the yml file in two places: First, in the `channels` section you will need to edit the path to the `ccdc_conda_channel` to reflect your own installation, and you will need to update the prefix at the bottom of the file to where you have installed the CSD Suite.

Once the yml has been updated, simply run:

```
conda env create -f STRIFE_environment.yml
conda activate STRIFE
```

## Setting Environment Variables

You will need to specify a handful of environment variables:

```
#Recommended - make these variables permanent by including in your .bashrc file
export CSDHOME=<path/to/CSD/installation>/CSD_2021 #Lets python see where the CSD installation is
export GHECOM_EXE=<path/to/ghecom/executable>/ghecom #Lets python know where the ghecom executable is for the buriedness calculations in the hotspots algorithm
export GOLD_DIR=<path/to/CSD/installation>/Discovery_2021/GOLD #To do GOLD docking in the Python API
```

After this step you should be able to run STRIFE from the terminal!


