#!/usr/bin/env python3.6

from main_module.modules import *


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
    parser.add_argument('-i', '--input', required=True,
                        help='Path to file with gene/features table')
    parser.add_argument('-ts', '--true_set', required=True,
                        help='Path to file with gene symbols of causal genes (True gene set)')
    parser.add_argument('-ass', '--algorithm_selection_set', required=False,
                        help='Path to file with algorithm selection set')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to output')
    parser.add_argument('--drop_qc', action='store_true', required=False,
                        help='(default=False) Perform pipeline without extended list and qc. \
                        Instead of optimal combination - use simple mean')
    parser.add_argument('--add_features', action='store_true', required=False,
                        help='(default=False) Add additional featutes to the provided table \
                        For detailed features description see github repo.')
    parser.add_argument('--set_seed', action='store_true', required=False,
                        help='(default=False) Switch it on if you want reproducibility')
    parser.add_argument('-n', '--n_bootstrap', type=int, required=False,
                        help='(default=15) Number of bootstraps for PU Bagging')
    parser.add_argument('-k', '--k_clusters', type=int, required=False,
                        help='(default=n_features/2) Number of clusters for Feature Agglomeration')
    
    args = parser.parse_args()
    
    input_file = import_fun.import_file(args.input)
    causal_genes = import_fun.import_file(args.true_set)

    
    if args.drop_qc == False:
        ass = import_fun.import_file(args.algorithm_selection_set)
        n_ass = import_fun.find_genes(ass.gene_symbol, input_file.gene_symbol)
        print(f'We found {n_ass} genes in ASS!')
    n_ts = import_fun.find_genes(causal_genes.gene_symbol, input_file.gene_symbol)
    print(f'We found {n_ts} genes from TS in provided file!')

    if args.add_features:
        input_file = proc_fts.add_features(input_file, causal_genes)

    X, y = proc_fun.return_x_y(input_file,
                               causal_genes,
                               n_clusters=import_fun.give_k(args.k_clusters, input_file))

    print('Launch ensemble...')

    ens = EnsembleClassifier(X, y, MODELS, ass, 
                             set_seed=args.set_seed,
                             n_bootstrap=import_fun.give_n(args.n_bootstrap))
    ens.run_estimators()
    probas = ens.best_scored_proba()
    probas.to_csv(f'{args.output}', sep='\t', index=False)
    print('Probabilities are returned...')
    print('Done!')
