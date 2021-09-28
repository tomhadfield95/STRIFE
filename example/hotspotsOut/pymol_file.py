
try:
    import tkinter as tk      
except ImportError:
    import Tkinter as tk
from os.path import join
import tempfile

import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()

dirpath = tempfile.mkdtemp()
zip_dir = "out.zip"
wd = os.getcwd()
with zipfile.ZipFile(zip_dir) as hs_zip:
    hs_zip.extractall(dirpath)

os.chdir(dirpath)
cmd.load("hotspot/apolar.ccp4", "apolar_hotspot")
cmd.isosurface(name="surface_apolar_hotspot", map="apolar_hotspot", level="5")

cmd.color("yellow", "surface_apolar_hotspot")
cmd.set("transparency", 0.2, "surface_apolar_hotspot")
cmd.load("hotspot/donor.ccp4", "donor_hotspot")
cmd.isosurface(name="surface_donor_hotspot", map="donor_hotspot", level="5")

cmd.color("blue", "surface_donor_hotspot")
cmd.set("transparency", 0.2, "surface_donor_hotspot")
cmd.load("hotspot/acceptor.ccp4", "acceptor_hotspot")
cmd.isosurface(name="surface_acceptor_hotspot", map="acceptor_hotspot", level="5")

cmd.color("red", "surface_acceptor_hotspot")
cmd.set("transparency", 0.2, "surface_acceptor_hotspot")
cmd.group("hotspot", members="apolar_hotspot")
cmd.group("hotspot", members="surface_apolar_hotspot")
cmd.group("hotspot", members="donor_hotspot")
cmd.group("hotspot", members="surface_donor_hotspot")
cmd.group("hotspot", members="acceptor_hotspot")
cmd.group("hotspot", members="surface_acceptor_hotspot")
cmd.load("hotspot/buriedness.ccp4", "buriedness_hotspot")
cmd.isosurface(name="surface_buriedness_hotspot", map="buriedness_hotspot", level="3")

cmd.color("gray", "surface_buriedness_hotspot")
cmd.set("transparency", 0.2, "surface_buriedness_hotspot")
cmd.group("hotspot", members="buriedness_hotspot")
cmd.group("hotspot", members="surface_buriedness_hotspot")
cmd.pseudoatom(object="PS_apolar_hotspot_0", pos=(5.0, 10.0, 4.0), color=(1, 1, 1), label=18.5)

cmd.pseudoatom(object="PS_apolar_hotspot_1", pos=(5.0, 10.0, 6.0), color=(1, 1, 1), label=19.3)

cmd.pseudoatom(object="PS_apolar_hotspot_2", pos=(6.100000000000001, 10.0, 2.5), color=(1, 1, 1), label=18.2)

cmd.pseudoatom(object="PS_apolar_hotspot_3", pos=(7.0, 10.0, 4.0), color=(1, 1, 1), label=19.3)

cmd.pseudoatom(object="PS_apolar_hotspot_4", pos=(-4.0, -8.5, -9.5), color=(1, 1, 1), label=17.8)

cmd.pseudoatom(object="PS_apolar_hotspot_5", pos=(-5.0, -10.5, -8.5), color=(1, 1, 1), label=18.4)

cmd.pseudoatom(object="PS_apolar_hotspot_6", pos=(-6.0, -10.0, -10.0), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_apolar_hotspot_7", pos=(10.5, 10.0, 1.0), color=(1, 1, 1), label=17.2)

cmd.pseudoatom(object="PS_apolar_hotspot_8", pos=(12.5, 12.0, 1.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_9", pos=(16.0, 7.5, -4.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_10", pos=(16.0, 9.0, -5.0), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_11", pos=(17.0, 7.0, -3.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_12", pos=(17.75, 9.5, -4.0), color=(1, 1, 1), label=15.8)

cmd.pseudoatom(object="PS_apolar_hotspot_13", pos=(18.0, 8.0, -3.0), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_14", pos=(19.0, 7.5, -2.0), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_apolar_hotspot_15", pos=(14.0, 21.5, -0.5), color=(1, 1, 1), label=15.1)

cmd.pseudoatom(object="PS_apolar_hotspot_16", pos=(15.0, 21.75, 2.5), color=(1, 1, 1), label=14.5)

cmd.pseudoatom(object="PS_apolar_hotspot_17", pos=(14.0, 22.0, 1.0), color=(1, 1, 1), label=15.1)

cmd.pseudoatom(object="PS_apolar_hotspot_18", pos=(13.5, 23.25, 3.25), color=(1, 1, 1), label=14.4)

cmd.pseudoatom(object="PS_apolar_hotspot_19", pos=(13.0, 23.0, 1.5), color=(1, 1, 1), label=15.1)

cmd.pseudoatom(object="PS_apolar_hotspot_20", pos=(14.5, 12.5, -1.0), color=(1, 1, 1), label=14.8)

cmd.pseudoatom(object="PS_apolar_hotspot_21", pos=(-7.5, -1.5, -14.0), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_apolar_hotspot_22", pos=(-6.5, -0.5, -10.0), color=(1, 1, 1), label=11.1)

cmd.pseudoatom(object="PS_apolar_hotspot_23", pos=(-6.5, -0.5, -12.0), color=(1, 1, 1), label=11.3)

cmd.pseudoatom(object="PS_apolar_hotspot_24", pos=(-3.0, -16.0, -8.0), color=(1, 1, 1), label=11.7)

cmd.pseudoatom(object="PS_apolar_hotspot_25", pos=(-1.0, -19.0, -6.5), color=(1, 1, 1), label=11.7)

cmd.pseudoatom(object="PS_apolar_hotspot_26", pos=(-1.5, -18.0, -7.5), color=(1, 1, 1), label=11.7)

cmd.pseudoatom(object="PS_apolar_hotspot_27", pos=(-2.25, -16.5, -9.75), color=(1, 1, 1), label=10.3)

cmd.pseudoatom(object="PS_apolar_hotspot_28", pos=(-2.5, -17.0, -7.0), color=(1, 1, 1), label=11.7)

cmd.pseudoatom(object="PS_apolar_hotspot_29", pos=(-1.0, -18.0, -9.0), color=(1, 1, 1), label=11.7)

cmd.pseudoatom(object="PS_apolar_hotspot_30", pos=(14.5, 7.0, 6.0), color=(1, 1, 1), label=11.2)

cmd.pseudoatom(object="PS_apolar_hotspot_31", pos=(15.5, 7.0, 3.5), color=(1, 1, 1), label=11.2)

cmd.pseudoatom(object="PS_apolar_hotspot_32", pos=(16.0, 7.0, 5.0), color=(1, 1, 1), label=11.2)

cmd.pseudoatom(object="PS_apolar_hotspot_33", pos=(17.0, 28.0, -0.5), color=(1, 1, 1), label=9.8)

cmd.pseudoatom(object="PS_apolar_hotspot_34", pos=(17.5, 28.0, 3.5), color=(1, 1, 1), label=9.8)

cmd.pseudoatom(object="PS_apolar_hotspot_35", pos=(17.0, 27.5, 2.0), color=(1, 1, 1), label=9.8)

cmd.pseudoatom(object="PS_apolar_hotspot_36", pos=(17.5, 28.5, 1.0), color=(1, 1, 1), label=9.8)

cmd.pseudoatom(object="PS_apolar_hotspot_37", pos=(19.0, 0.0, 4.0), color=(1, 1, 1), label=8.4)

cmd.pseudoatom(object="PS_apolar_hotspot_38", pos=(20.0, -2.5, 7.5), color=(1, 1, 1), label=8.4)

cmd.pseudoatom(object="PS_apolar_hotspot_39", pos=(20.0, -0.75, 6.0), color=(1, 1, 1), label=6.9)

cmd.pseudoatom(object="PS_apolar_hotspot_40", pos=(19.25, -2.75, 5.5), color=(1, 1, 1), label=7.2)

cmd.pseudoatom(object="PS_apolar_hotspot_41", pos=(21.5, 8.5, -8.0), color=(1, 1, 1), label=8.2)

cmd.pseudoatom(object="PS_apolar_hotspot_42", pos=(17.833333333333336, -2.25, -5.666666666666664), color=(1, 1, 1), label=4.8)

cmd.pseudoatom(object="PS_apolar_hotspot_43", pos=(18.0, -3.0, -7.5), color=(1, 1, 1), label=7.5)

cmd.pseudoatom(object="PS_apolar_hotspot_44", pos=(16.5, -2.5, -3.5), color=(1, 1, 1), label=7.5)

cmd.pseudoatom(object="PS_apolar_hotspot_45", pos=(16.25, -3.25, -6.5), color=(1, 1, 1), label=4.3)

cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_0")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_0")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_1")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_1")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_2")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_2")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_3")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_3")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_4")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_4")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_5")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_5")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_6")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_6")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_7")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_7")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_8")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_8")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_9")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_9")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_10")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_10")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_11")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_11")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_12")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_12")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_13")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_13")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_14")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_14")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_15")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_15")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_16")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_16")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_17")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_17")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_18")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_18")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_19")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_19")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_20")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_20")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_21")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_21")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_22")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_22")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_23")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_23")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_24")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_24")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_25")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_25")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_26")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_26")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_27")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_27")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_28")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_28")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_29")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_29")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_30")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_30")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_31")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_31")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_32")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_32")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_33")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_33")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_34")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_34")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_35")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_35")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_36")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_36")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_37")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_37")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_38")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_38")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_39")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_39")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_40")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_40")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_41")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_41")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_42")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_42")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_43")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_43")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_44")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_44")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_45")
cmd.group("label_apolar_hotspot", members="PS_apolar_hotspot_45")
cmd.pseudoatom(object="PS_donor_hotspot_0", pos=(-1.0, -19.5, -10.0), color=(1, 1, 1), label=30.1)

cmd.pseudoatom(object="PS_donor_hotspot_1", pos=(14.5, 9.5, -4.5), color=(1, 1, 1), label=19.3)

cmd.pseudoatom(object="PS_donor_hotspot_2", pos=(15.0, 26.0, -1.5), color=(1, 1, 1), label=19.1)

cmd.pseudoatom(object="PS_donor_hotspot_3", pos=(3.5, 10.0, 1.5), color=(1, 1, 1), label=19.0)

cmd.pseudoatom(object="PS_donor_hotspot_4", pos=(8.0, 9.0, 6.0), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_donor_hotspot_5", pos=(-2.0, -9.0, -9.0), color=(1, 1, 1), label=18.0)

cmd.pseudoatom(object="PS_donor_hotspot_6", pos=(14.5, 6.5, 0.5), color=(1, 1, 1), label=17.3)

cmd.pseudoatom(object="PS_donor_hotspot_7", pos=(19.0, 10.0, -3.0), color=(1, 1, 1), label=16.5)

cmd.pseudoatom(object="PS_donor_hotspot_8", pos=(15.0, 6.5, -4.5), color=(1, 1, 1), label=16.4)

cmd.pseudoatom(object="PS_donor_hotspot_9", pos=(13.0, 7.5, 1.5), color=(1, 1, 1), label=16.1)

cmd.pseudoatom(object="PS_donor_hotspot_10", pos=(-4.5, 0.0, -12.5), color=(1, 1, 1), label=16.0)

cmd.pseudoatom(object="PS_donor_hotspot_11", pos=(17.0, 10.0, -7.0), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_donor_hotspot_12", pos=(15.0, 9.5, -2.0), color=(1, 1, 1), label=15.6)

cmd.pseudoatom(object="PS_donor_hotspot_13", pos=(10.5, 7.5, 3.5), color=(1, 1, 1), label=15.5)

cmd.pseudoatom(object="PS_donor_hotspot_14", pos=(16.0, 25.0, 3.0), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_donor_hotspot_15", pos=(9.5, 8.0, 1.0), color=(1, 1, 1), label=14.1)

cmd.pseudoatom(object="PS_donor_hotspot_16", pos=(14.0, 9.0, 4.5), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_donor_hotspot_17", pos=(-1.0, -19.0, -7.0), color=(1, 1, 1), label=13.2)

cmd.pseudoatom(object="PS_donor_hotspot_18", pos=(17.0, 24.0, 0.0), color=(1, 1, 1), label=13.0)

cmd.pseudoatom(object="PS_donor_hotspot_19", pos=(14.0, 11.0, -6.0), color=(1, 1, 1), label=12.7)

cmd.pseudoatom(object="PS_donor_hotspot_20", pos=(12.5, 9.5, 6.5), color=(1, 1, 1), label=11.0)

cmd.pseudoatom(object="PS_donor_hotspot_21", pos=(17.0, 6.0, -1.0), color=(1, 1, 1), label=10.7)

cmd.pseudoatom(object="PS_donor_hotspot_22", pos=(17.5, -1.5, 5.5), color=(1, 1, 1), label=10.6)

cmd.pseudoatom(object="PS_donor_hotspot_23", pos=(-2.5, 12.0, -13.0), color=(1, 1, 1), label=9.4)

cmd.pseudoatom(object="PS_donor_hotspot_24", pos=(19.0, 0.0, 4.0), color=(1, 1, 1), label=9.3)

cmd.pseudoatom(object="PS_donor_hotspot_25", pos=(0.0, 16.0, -11.0), color=(1, 1, 1), label=9.2)

cmd.pseudoatom(object="PS_donor_hotspot_26", pos=(17.0, 28.5, -4.5), color=(1, 1, 1), label=8.5)

cmd.pseudoatom(object="PS_donor_hotspot_27", pos=(14.0, 8.0, 11.0), color=(1, 1, 1), label=7.8)

cmd.pseudoatom(object="PS_donor_hotspot_28", pos=(-3.0, -16.5, -7.0), color=(1, 1, 1), label=7.8)

cmd.pseudoatom(object="PS_donor_hotspot_29", pos=(-3.0, 11.0, -9.5), color=(1, 1, 1), label=7.8)

cmd.pseudoatom(object="PS_donor_hotspot_30", pos=(19.0, 26.5, -3.0), color=(1, 1, 1), label=6.8)

cmd.pseudoatom(object="PS_donor_hotspot_31", pos=(18.5, 22.0, 3.0), color=(1, 1, 1), label=6.4)

cmd.group("label_donor_hotspot", members="PS_donor_hotspot_0")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_0")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_1")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_1")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_2")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_2")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_3")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_3")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_4")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_4")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_5")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_5")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_6")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_6")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_7")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_7")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_8")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_8")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_9")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_9")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_10")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_10")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_11")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_11")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_12")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_12")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_13")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_13")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_14")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_14")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_15")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_15")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_16")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_16")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_17")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_17")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_18")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_18")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_19")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_19")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_20")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_20")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_21")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_21")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_22")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_22")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_23")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_23")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_24")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_24")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_25")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_25")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_26")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_26")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_27")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_27")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_28")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_28")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_29")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_29")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_30")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_30")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_31")
cmd.group("label_donor_hotspot", members="PS_donor_hotspot_31")
cmd.pseudoatom(object="PS_acceptor_hotspot_0", pos=(3.0, 10.5, 4.5), color=(1, 1, 1), label=21.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_1", pos=(18.0, 11.5, -2.5), color=(1, 1, 1), label=18.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_2", pos=(8.0, 10.0, 2.5), color=(1, 1, 1), label=18.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_3", pos=(-3.5, -7.0, -9.5), color=(1, 1, 1), label=18.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_4", pos=(-4.0, -9.0, -11.0), color=(1, 1, 1), label=18.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_5", pos=(8.0, 8.5, 5.0), color=(1, 1, 1), label=17.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_6", pos=(7.5, 8.0, 2.5), color=(1, 1, 1), label=17.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_7", pos=(16.0, 9.5, -5.5), color=(1, 1, 1), label=15.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_8", pos=(9.0, 8.5, 1.0), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_9", pos=(14.0, 25.5, 1.0), color=(1, 1, 1), label=15.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_10", pos=(16.0, 20.0, 3.0), color=(1, 1, 1), label=15.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_11", pos=(18.5, 5.0, -2.0), color=(1, 1, 1), label=14.3)

cmd.pseudoatom(object="PS_acceptor_hotspot_12", pos=(-5.5, -2.0, -9.0), color=(1, 1, 1), label=13.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_13", pos=(14.0, 11.5, -4.5), color=(1, 1, 1), label=13.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_14", pos=(15.5, 13.5, -1.5), color=(1, 1, 1), label=12.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_15", pos=(-9.5, -11.5, -10.5), color=(1, 1, 1), label=12.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_16", pos=(12.5, 13.5, 1.5), color=(1, 1, 1), label=12.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_17", pos=(-9.5, -10.0, -12.5), color=(1, 1, 1), label=12.2)

cmd.pseudoatom(object="PS_acceptor_hotspot_18", pos=(-2.5, -15.0, -10.0), color=(1, 1, 1), label=11.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_19", pos=(18.5, -1.5, -3.5), color=(1, 1, 1), label=10.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_20", pos=(18.0, 22.5, 2.0), color=(1, 1, 1), label=10.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_21", pos=(-0.5, 16.5, -11.0), color=(1, 1, 1), label=10.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_22", pos=(16.5, 28.5, -3.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_23", pos=(0.0, -15.0, -14.5), color=(1, 1, 1), label=10.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_24", pos=(-0.5, -17.0, -10.0), color=(1, 1, 1), label=9.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_25", pos=(16.5, 6.0, 2.5), color=(1, 1, 1), label=9.5)

cmd.pseudoatom(object="PS_acceptor_hotspot_26", pos=(1.0, -17.5, -9.0), color=(1, 1, 1), label=9.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_27", pos=(20.5, 10.0, -6.5), color=(1, 1, 1), label=9.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_28", pos=(22.0, 10.5, -2.5), color=(1, 1, 1), label=8.8)

cmd.pseudoatom(object="PS_acceptor_hotspot_29", pos=(-3.5, 0.5, -13.5), color=(1, 1, 1), label=8.1)

cmd.pseudoatom(object="PS_acceptor_hotspot_30", pos=(15.5, -3.5, -7.0), color=(1, 1, 1), label=7.9)

cmd.pseudoatom(object="PS_acceptor_hotspot_31", pos=(7.0, 9.0, -1.5), color=(1, 1, 1), label=7.7)

cmd.pseudoatom(object="PS_acceptor_hotspot_32", pos=(14.0, 10.0, 6.5), color=(1, 1, 1), label=7.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_33", pos=(22.0, 8.0, -9.5), color=(1, 1, 1), label=6.6)

cmd.pseudoatom(object="PS_acceptor_hotspot_34", pos=(19.0, -4.0, -7.0), color=(1, 1, 1), label=6.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_35", pos=(9.0, 10.5, -2.5), color=(1, 1, 1), label=6.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_36", pos=(13.0, 10.0, -2.0), color=(1, 1, 1), label=6.4)

cmd.pseudoatom(object="PS_acceptor_hotspot_37", pos=(5.0, -12.5, -14.0), color=(1, 1, 1), label=6.0)

cmd.pseudoatom(object="PS_acceptor_hotspot_38", pos=(11.5, 8.5, 10.0), color=(1, 1, 1), label=5.3)

cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_0")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_0")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_1")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_1")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_2")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_2")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_3")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_3")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_4")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_4")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_5")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_5")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_6")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_6")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_7")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_7")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_8")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_8")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_9")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_9")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_10")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_10")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_11")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_11")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_12")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_12")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_13")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_13")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_14")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_14")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_15")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_15")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_16")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_16")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_17")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_17")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_18")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_18")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_19")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_19")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_20")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_20")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_21")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_21")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_22")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_22")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_23")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_23")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_24")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_24")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_25")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_25")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_26")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_26")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_27")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_27")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_28")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_28")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_29")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_29")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_30")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_30")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_31")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_31")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_32")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_32")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_33")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_33")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_34")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_34")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_35")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_35")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_36")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_36")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_37")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_37")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_38")
cmd.group("label_acceptor_hotspot", members="PS_acceptor_hotspot_38")
cmd.group("labels_hotspot", members="label_apolar_hotspot")
cmd.group("labels_hotspot", members="label_donor_hotspot")
cmd.group("labels_hotspot", members="label_acceptor_hotspot")
cmd.load("hotspot/protein.pdb", "protein_hotspot")
cmd.show("cartoon", "protein_hotspot")
cmd.hide("line", "protein_hotspot")
cmd.show("sticks", "organic")


class IsoLevel(tk.Variable):
    def __init__(self, master, name, level):
        tk.Variable.__init__(self, master, value=level)
        self.name = name
        self.trace('w', self.callback)

    def callback(self, *args):
        cmd.isolevel(self.name, self.get())

    def increment(self, event=None, delta=0.1):
        self.set(round(float(self.get()) + delta, 2))

    def decrement(self, event=None):
        self.increment(None, -0.1)


surface_list = {'hotspot': {'fhm': ['surface_apolar_hotspot', 'surface_donor_hotspot', 'surface_acceptor_hotspot'], 'buriedness': ['surface_buriedness_hotspot']}}
surface_max_list = {'hotspot': {'fhm': 30.1, 'buriedness': 8}}

top = tk.Toplevel(plugins.get_tk_root())

master = tk.Frame(top, padx=10, pady=10)
master.pack(fill="both", expand=1)

for child in list(master.children.values()):
    child.destroy()


row_counter = 0
for identifier, component_dic in surface_list.items():
    # add calculation identifier
    tk.Label(master, text=identifier).grid(row=row_counter, column=0, sticky="w")
    row_counter += 1
    
    for component_id, surfaces in component_dic.items():
        # add collection label, e.g. superstar or hotspot etc.
        tk.Label(master, text=component_id).grid(row=row_counter, column=1, sticky='w')
        row_counter += 1
        
        for i, surface in enumerate(surfaces):
            # add grid type label
            probe = surface.split("_")[-2]
            tk.Label(master, text=probe).grid(row=row_counter, column=2, sticky="w")
            
            # slider code 
            v = IsoLevel(master, surface, 5)
            e = tk.Scale(master, orient=tk.HORIZONTAL, from_=0, to=surface_max_list[identifier][component_id],
                         resolution=0.1, showvalue=0, variable=v)
            e.grid(row=row_counter, column=3, sticky="ew")

            e = tk.Entry(master, textvariable=v, width=4)
            e.grid(row=row_counter, column=4, sticky="e")
            master.columnconfigure(3, weight=1)
            row_counter += 1



cmd.bg_color("white")
if wd:
    os.chdir(wd)