import os.path
import pandas as pd
def loadEgPed():
    '''Load example pedigree data.

    '''
    basepath = os.path.abspath(__file__)
    folder = os.path.dirname(basepath)
    data_path = os.path.join(folder, 'data/ped.txt')
    text = pd.read_table(data_path,header=0)

    return text
def loadEgGeno():
    '''Load example genotype data.

    '''
    basepath = os.path.abspath(__file__)
    folder = os.path.dirname(basepath)
    data_path = os.path.join(folder, 'data/geno.traw')


    return data_path 
