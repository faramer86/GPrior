# Gene Prioritization tool (GPrior)

GPrior is a Gene Prioritization tool that uses PU learning to prioritize candidate genes based on their similarity with a set of known positive examples. As a result it returns a table with probabilities for each gene.

In order to perform prioritization user can repeat our pipeline from the article or compile his own table with features.

## Dependencies

Make sure that you have `python` version >= 3.6 and latest version of `pip`. You can check it by running:

```bash
 $ python --version
 $ pip --version
```

If you do not have Python, please install the latest 3.x version from python.org.

In order GPrior to be cross-platform we use `virtualenv`. If you do not have it you can install it using `pip` (use `sudo` if necessary):

```bash
pip install virtualenv
```

If you have any problems with the installation or usage  of `virtualenv` be free to consult the official documentation: https://virtualenv.pypa.io/en/latest/

Go to the GPrior repository, create and launch virtual environment:

```
virtualenv vprior
source vprior/bin/activate
```

If you see `(vprior)` prefix in your terminal, then everything goes well so far.

Install the requirments from the home repository:

```
pip install -r requirements.txt
```

This command will install all the necessary python packages. You can go to the launching `.gprior.py`, it should work.

Do not forget to exit virtual environment after work:

```
deactivate
```

## Usage

#### gprior.py

```bash
  -h, --help            show this help message and exit
  -i, --input
                        Path to file with gene/features table
  -ts, --true_set 
                        Path to file with gene symbols of causal genes (True
                        gene set)
  -ass, --algorithm_selection_set 
                        Path to file with algorithm selection set
  -o, --output
                        Path to output
  --drop_qc             Perform pipeline without extended list and qc. Instead
                        of optimal combination - use simple mean
  --add_features        Add additional featutes to the provided table For
                        detailed features description see github repo.
  --set_seed            Switch it on if you want reproducibility

```

example:

```bash
 ./gprior.py \ # file that you launch
    -i test_input_copy.tsv \ # input .tsv file with genes and features (see input example) 
    -ts test_ts.tsv \ # True set of genes (one column - "gene_symbol")
    -ass test_ass.tsv \ # Algorithm selection set of genes (one column - "gene_symbol")
    -o test_output.tsv # output
```



## Arguments


## License
[MIT](https://choosealicense.com/licenses/mit/)
