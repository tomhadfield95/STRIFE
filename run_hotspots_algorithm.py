#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 16:56:34 2021

@author: hadfield
"""

from preprocessing import preprocessing
import argparse


def hotspots_calculation(p, args):
    
    p.calculateFHM(args.pdb, args.output_dir)

if __name__=='__main__':
        
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str, help='Location of Protein pdb')
    parser.add_argument('output_dir', type=str,
                        help='Directory in which to store outputs')
    
    arguments = parser.parse_args()
    
    p = preprocessing()
    
    hotspots_calculation(p, arguments)
    
    