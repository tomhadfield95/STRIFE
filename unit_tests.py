#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:37:28 2022

@author: hadfield

Script which tests whether STRIFE has successfully installed 
    Checks that csd_python_api and hotspots are installed in the conda environment
    And that key environment variables have been set

"""

import unittest
import os
import glob


class TestSTRIFE(unittest.TestCase):
    
    
    def test_csdhome(self):
        print('\n')
        print('Checking that CSDHOME shell variable has been correctly set')
    
        CSDHOME =os.getenv('CSDHOME')
        
        self.assertFalse(CSDHOME == '', msg='CSDHOME variable has not been set')
        self.assertTrue(os.path.exists(CSDHOME), 'CSDHOME variable points to a directory that does not exist')
        self.assertTrue(os.path.exists(os.path.join(CSDHOME, 'bin', 'activate')), 'CSDHOME variable points to a directory that does not contain an installation of the CSD API')
    
    
    def test_gold_dir(self):
        
        print('\n')
        print('Checking that GOLD_DIR shell variable has been correctly set')
        
        GOLD_DIR = os.getenv('GOLD_DIR')
        
        self.assertFalse(GOLD_DIR == '', msg='GOLD_DIR variable has not been set')
        self.assertTrue(os.path.exists(GOLD_DIR), 'GOLD_DIR variable points to a directory that does not exist')
        self.assertTrue(os.path.exists(os.path.join(GOLD_DIR, 'bin', 'gold_auto')), 'GOLD_DIR variable points to a directory that does not contain an installation of GOLD')
        
    def test_ghecom_exe(self):
        
        print('\n')
        print('Checking that GHECOM_EXE shell variable has been correctly set')
        
        GHECOM_EXE = os.getenv('GHECOM_EXE')
        
        self.assertFalse(GHECOM_EXE == '', msg='GHECOM_EXE variable has not been set')
        self.assertTrue(os.path.exists(GHECOM_EXE), 'GHECOM_EXE variable does not point to an existing file')
        self.assertTrue(GHECOM_EXE[-6:] == 'ghecom', 'GHECOM_EXE does not point to an executable file named ghecom')
        
def kw_in_list(list_of_strings, kw):
    
    present = False
    
    for string in list_of_strings:
        if kw in string:
            present = True
    
    return present
        
if __name__=='__main__':
    print('Attempting to import important modules')
    

    all_packages = glob.glob(f'{os.getenv("CONDA_PREFIX")}/lib/python3.7/site-packages/*')
    all_packages_suffix = [x.split('/')[-1] for x in all_packages]

    csd_api_installed = kw_in_list(all_packages_suffix, 'csd_python_api')
    hotspots_api_installed = kw_in_list(all_packages_suffix, 'hotspots')

    
    if not csd_api_installed:
        raise ModuleNotFoundError('CSD API does not appear to be installed in this conda environment!')
    
    if not hotspots_api_installed:
        raise ModuleNotFoundError('Hotspots API does not appear to be installed in this conda environment! Please check the Github README for guidance on how to install it')

    print('Now run unit test')
    
    unittest.main()

    