#!/usr/bin/env python3.6

from modules import *
from density_plots import create_density_plot

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
    parser.add_argument('-w', '--class_weights',
                        required=False, default={1: 1, 0: 100},
                        help='Class weights in model')
    parser.add_argument('-m', '--model',
                        required=False, default='rf',
                        help='Choose appropriate algorithm for classification')
    args = parser.parse_args()
    input = import_postgap_file(args.postgap_path, COL_DTYPES)
    causal_genes = pd.read_csv(args.genes, sep='\t')
    print('List of causal genes is imported...')
    print('Processing input...')
    ML_input = process_input_file(input, causal_genes, GTEX_COLUMNS, GTEX_DB)
    print('Processing input is done...')
    # print('Launch classifier...')
    # probabilities = GeneClassifier(
    #     ML_input, causal_genes, args.model, args.class_weights) \
    #     .return_probability()
    # print('Probabilities are returned...')
    # print(f'Writing probabilities for genes to be causal to "{args.output}"')
    ML_input.to_csv(Path(args.output), sep='\t', index=False)
    # print('Create density plots of gene names from extended list...')
    # if args.extended_list:
    #     extended_list = pd.read_csv(args.extended_list, sep='\t')
    #     print('Extended list of genes is imported...')
    #     create_density_plot(probabilities, extended_list,
    #                         os.path.split(args.output)[0])
    print('Done!')
