#!/usr/bin/env python3.6

from modules import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='This script analyse gwas data')
    parser.add_argument('-i', '--postgap_path', required=True,
                        help='Path to file/folder with postgap output')
    parser.add_argument('-g', '--genes', required=True,
                        help='Path to file with names of causal genes')
    parser.add_argument('-e', '--extended_list', required=True,
                        help='Path to file with extended list of causal genes')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to output')
    args = parser.parse_args()

    input = import_postgap_file(args.postgap_path)
    causal_genes = pd.read_csv(args.genes, sep='\t')
    print('List of causal genes is imported...')
    extended_list = pd.read_csv(args.extended_list, sep='\t')
    print('Extended list of genes is imported...')
    print('Processing input...')
    ML_input = process_input_file(input, gtex_columns, gtex_db)
    print('Processing input is done...')
    # ML part:
    print('Launch classifier...')
    probabilities = GeneClassifier(
        input_df, causal_genes, extended_list).return_probability()
    print('Probabilities are returned...')
    print(f'Write probabilities for genes to be causal to "{args.output}"')
    ML_input.to_csv(Path(args.output), sep='\t', index=False)
    print('Done!')
