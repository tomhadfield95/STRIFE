#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 12:25:12 2021

@author: hadfield
"""

#SCRIPT MUST BE RUN IN A CONDA ENVIRONMENT WHICH SUPPORTS PYMOL

from pymol import cmd, finish_launching, plugins
from pymol.cgo import *
import sys

outStoreDir = sys.argv[1]

finish_launching()



#Load into pymol
cmd.load(f'{outStoreDir}/protein.pdb', "prot") #Load protein
cmd.show("cartoon", "prot")
cmd.hide("line", "prot")
cmd.show("sticks", "organic")

cmd.load(f'{outStoreDir}/frag.sdf', 'fragment') #Load fragment
cmd.load(f'{outStoreDir}/pharmacophoreLattice.sdf', 'pharm_lattice') #Load lattice

cmd.set('all_states', 'on') #Change pymol settings so that all atoms in the lattice can be viewed at once.
#cmd.set("sphere_scale", 0.035, "pharm_lattice") #Change size of lattice size atoms for easier viewability
cmd.set("sphere_scale", 0.05, "pharm_lattice") #Change size of lattice size atoms for easier viewability
cmd.set("bg_rgb", [1,1,1])

cmd.cd(dir = f'{outStoreDir}')

