#!/usr/bin/env python3.6

from gprior.modules import *


warnings.filterwarnings("ignore", category=FutureWarning)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='./gprior.py',
                            usage='%(prog)s [options]',
                            description="""
                            GPrior is a Gene Prioritization tool that uses PU learning
                            to prioritize candidate genes based on their similarity with
                            a set of known positive examples. As a result it returns a
                            table with probabilities for each gene.
                            """)
    parser.add_argument('-i', '--input', required=True,  metavar='',
                        help='Path to gene-features table')
    parser.add_argument('-ts', '--true_set', required=True, metavar='',
                        help='Path to table with gene symbols of causal genes (True gene set)')
    parser.add_argument('-aes', '--algorithm_evaluation_set', required=False,metavar='',
                        help='Path to file with algorithm evaluation set')
    parser.add_argument('-o', '--output', required=True, metavar='',
                        help='Path to output file')
    parser.add_argument('-n', '--n_bootstrap', type=int, required=False, metavar='',
                        help='(default=15) Number of bootstraps for PU Bagging')
    parser.add_argument('-k', '--k_clusters', type=int, required=False, metavar='',
                        help='(default=n_features) Number of clusters for Feature Agglomeration')
    parser.add_argument('-s', '--s_coef', type=float, required=False, metavar='',
                        help='(default=1) Sampling coefficient')
    parser.add_argument('--drop_aes', action='store_true', required=False,
                        help='(positional; default=False) Perform pipeline without AES and qc. \
                        Instead of optimal combination - simple mean will be used')
    parser.add_argument('--add_features', action='store_true', required=False,
                        help='(positional; default=False) Add additional featutes to the provided table')
    parser.add_argument('--tune', action='store_true', required=False,
                        help='(positional; default=False) Tune hyperparameters of the model \
                        It significently slow down training process. Do not use it with high n.')
    parser.add_argument('--set_seed', action='store_true', required=False,
                        help='(positional; default=False) Switch it on if you want reproducibility')                  
    
    args = parser.parse_args()
    
    input_file = INtoolbox.import_file(args.input)
    causal_genes = INtoolbox.import_file(args.true_set)

    if args.algorithm_evaluation_set:
        aes = INtoolbox.import_file(args.algorithm_evaluation_set)
        n_aes = INtoolbox.find_genes(aes.gene_symbol, input_file.gene_symbol)
        print(f'There are {n_aes} genes from AES found in input file!')
    else:
        print('NOTE: You have not specified AES set!',
              'NOTE: Average prediction will be calculated!', sep='\n')
    
    n_ts = INtoolbox.find_genes(causal_genes.gene_symbol, input_file.gene_symbol)
    print(f'There are {n_ts} genes from TS found in input file!')

    if args.add_features:
        input_file = proc_fts.add_features(input_file, causal_genes)

    X, y = MLtoolbox.return_x_y(input_file,
                               causal_genes,
                               k_clusters=INtoolbox.give_k(args.k_clusters, input_file))

    ens = EnsembleClassifier(X, y, MODELS, 
                             set_seed=args.set_seed,
                             tune=args.tune,
                             n_bootstrap=INtoolbox.give_n(args.n_bootstrap),
                             s_coef=INtoolbox.give_s(args.s_coef))
    
    if not args.drop_aes and args.algorithm_evaluation_set:
        ens.set_aes(aes)
        ens.set_ytrue()
        ens.run_estimators()
        probas = ens.best_scored_proba()
    else:
        ens.run_estimators()
        probas = ens.get_prediction_set()

    probas.to_csv(args.output, sep='\t', index=False)

    print('Done!')
