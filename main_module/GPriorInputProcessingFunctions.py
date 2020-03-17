from main_module.var import *
from main_module.AdditionalFeatures import *
import pandas as pd
import os
import glob

def find_genes(list_1, list_2):
    return len(set(list_1) & set(list_2))


def import_file(input_path):
    """
    Import files (--input, -ts, -ass)
    """
    filename = os.path.basename(input_path)

    assert '.tsv' in filename, 'Plese, use .tsv extansion for all of the files!'
    
    if os.path.isfile(input_path):
        print(f'Importing {filename}')
        infile = pd.read_csv(input_path, sep='\t', low_memory=False)
        if 'gene_symbol' not in list(infile):
            sys.exit(f'({filename}) Please specify column with genes in all of the files as "gene_symbol"!')
        return infile
    else:
        sys.exit(f'There is no file with such name ({filename}) or it is not a file!')
