import os, sys
import numpy as np
import pandas as pd

def reverse_complement(seq):
    comp_dict = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'N': 'N'
    }
    seq_list = list(seq.upper()[::-1])
    comp_arr = [comp_dict[seq_list[0]],comp_dict[seq_list[1]],comp_dict[seq_list[2]],comp_dict[seq_list[3]]]
    
    return ''.join(comp_arr)
