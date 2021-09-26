#To run the STRIFE algorithm you need to export the following variables:

#Recommended - make these variables permanent by including in your .bashrc file
export CSDHOME=<path/to/CSD/installation>/CSD_2021 #Lets python see where the CSD installation is
export GHECOM_EXE=<path/to/ghecom/executable>/ghecom #Lets python know where the ghecom executable is for the buriedness calculations in the hotspots algorithm
export GOLD_DIR=<path/to/CSD/installation>/Discovery_2021/GOLD #To do GOLD docking in the Python API



#WARNING - your LD_LIBRARY_PATH should usually look like it is defined in the line below. However, it needs to be changed when running the STRIFE algorithm. Unfortunately, on some machines when you change LD_LIBRARY_PATH to the definition two lines below, it seems to break the 'ssh' functionality on your computer. This is obviously very annoying, but it can be easily gotten around by opening a new terminal tab/window or redefining your LD_LIBARY_PATH to look like the firsone below 

export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib:/usr/lib/X11 #Should normally look like this (or similar)
export LD_LIBRARY_PATH=<path/to/CSD/installation>/Python_API_2021/miniconda/lib:<path/to/CSD/installation>/Python_API_2021/miniconda/lib/python3.7/site-packages/ccdc/_lib/:/lib/:/usr/lib:/usr/local/lib:/usr/lib/X11 #Should be changed to look like this before running the STRIFE algorithm
