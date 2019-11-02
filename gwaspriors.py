#!/usr/bin/env python3.6

from modules import *
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='./gwaspriors.py',
                            usage='%(prog)s [options]',
                            description="""
                            This script prioritize genes
                            using gwas summary statistics.
                            It has two modes:
                            First, data preprocessing (use "--merge") 
                            (Merging postgap data;
                            Adding new features).
                            Second, gene prioritization
                            (PU-learning using ensemble classifier).
                            """)
    parser.add_argument('-i', '--input', required=False,
                        help='Path to file/folder with postgap output')
    parser.add_argument('-g', '--genes', required=True,
                        help='Path to file with names of causal genes')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to output')
    parser.add_argument('-e', '--extended_list', required=False,
                        help='Path to file with extended list of causal genes')
    parser.add_argument('--merge', action='store_true',
                        help='Switch to preprocessing mode')
    args = parser.parse_args()

    causal_genes = pd.read_csv(args.genes, sep='\t')

    if args.merge:
        input_file = import_postgap_file(args.input, COL_DTYPES)
        process_input_file(input_file, causal_genes).to_csv(Path(args.output), sep='\t')
        sys.exit("Done!")
    
    input_file = pd.read_csv(args.input, sep='\t')
    extended_list = pd.read_csv(args.extended_list, sep='\t')
    n_genes_in_extended = number_of_intersections(extended_list['gene_symbol'], 
                                                    input_file['gene_symbol'])
    print(f'We found {n_genes_in_extended} genes in extended list!')
    X, y = return_x_y(input_file, causal_genes)
    print('Launch ensemble...\n')
    ens = EnsembleClassifier(X, y, MODELS, extended_list)
    ens.run_estimators()
    probas = ens.best_scored_proba()
    probas.to_csv(f'{args.output}', sep='\t', index=False)
    print('Probabilities are returned...')
    print('Done!')
