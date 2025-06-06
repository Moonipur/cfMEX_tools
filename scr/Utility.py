"""
Code: Utility (Version 1.2)
By: Songphon Sutthitthasakul (Moon)
Date: 06-06-2025
"""

import os, sys
import numpy as np
import pandas as pd

def reverse_complement(seq):
    comp_dict = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N',
    'R': 'N',
    'Y': 'N',
    'K': 'N',
    'M': 'N',
    'W': 'N',
    'S': 'N',
    'B': 'N',
    'D': 'N',
    'H': 'N',
    'V': 'N'
    }
    seq_list = list(seq.upper()[::-1])
    comp_arr = [comp_dict[seq_list[0]],comp_dict[seq_list[1]],comp_dict[seq_list[2]],comp_dict[seq_list[3]]]
    
    return ''.join(comp_arr)
