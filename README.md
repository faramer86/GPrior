![](logo.png)

# Gene Prioritization tool (GPrior)

GPrior is a Gene Prioritization tool that uses PU learning to prioritize candidate genes based on their similarity with a set of known positive examples. As a result it returns a table with probabilities for each gene.

## Table of Contents

1. [Usage](#usage)
2. [Dependencies](#dependencies)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Features](#features)
6. [Input example](#input-example)
7. [Output example](#output-example)

## Usage

In order to perform prioritization user can repeat our pipeline from the article or compile his own table with features.

## Dependencies 

Make sure that you have `python` version >= 3.6 and latest version of `pip`. You can check it by running:

```bash
 python --version
 pip --version
```

If you do not have Python, please install the latest 3.x version from python.org.

In order GPrior to be cross-platform we use `virtualenv`. If you do not have it you can install it using `pip` (use `sudo` if necessary):

```bash
pip install virtualenv
```

If you have any problems with the installation or usage  of `virtualenv` be free to consult the official documentation: https://virtualenv.pypa.io/en/latest/

## Installation

Firstly, clone GPrior repository to your local machine. 

```bash
git clone https://github.com/faramer86/GPrior.git
```

Go to the GPrior repository, create and launch virtual environment:

```bash
virtualenv vprior
source vprior/bin/activate
```

If you see `(vprior)` prefix in your terminal, then everything goes well so far.

Install the requirments from the home repository. 
This command will install all the necessary python packages:

```bash
pip install -r requirements.txt
```

Now you can launch `gprior.py`, it should work.

#### NOTE:

1) Do not forget to exit virtual environment after you are done:

```bash
deactivate
```

2) Do not forget to activate `vprior` everytime you want to launch tool:

```bash
source vprior/bin/activate
```

## Usage

There are two main scripts: `gprior.py` and `process_postgap.py`. If you want to prioritize genes based on your own table of features - use `gprior.py`. If you want to reproduce pipeline from the article - use `process_postgap.py`. For more detailes read futher.

### gprior.py

`gprior.py` uses table of features and sets of causal genes as an input (see **Input examples** section) and perform gene prioritization. It consists of 5 PU Bagging classifiers. All the predictions sre combined using optimal combination approach (see article). True gene set (TS) is used for training and algorithm selection set (ASS) for optimal combination and quality evaluation.

We add several useful tunable functional features. See `--help` page:

### Arguments:

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

### process_postgap.py

## Features

## Input example:

1) `./gprior.py` Input file (`-i`):

 |**gene_symbol**| **feature 1** |**feature 2**|**...**|**feature n**|
 |:----:| :--------------------: |:--------------------:|---|:--------------------:|
 |GENE SYMBOL 1| ... </br>  | ... </br> |...|... </br> |
 |...|...| ... |...|...|...|
 |GENE SYMBOL n| ... </br> |... </br> |...|... </br>|

2) `./process_postgap.py` Input file (`-i`):

 |**gene_symbol**| ** SNPs ** |**feature 2**|**...**|**feature n**|
 |:----:| :--------------------: |:--------------------:|---|:--------------------:|
 |GENE SYMBOL 1| SNP1 </br>  | ... </br> |...|... </br> |
 |GENE SYMBOL 1| SNP2 </br>  | ... </br> |...|... </br> |
 |GENE SYMBOL 1| SNP3 </br>  | ... </br> |...|... </br> |
 |...|...| ... |...|...|...|
 |GENE SYMBOL n| ... </br> |... </br> |...|... </br>|

For more detailes see: https://github.com/Ensembl/postgap

3) Algorithm selection set (`-ass`) and True gene set (`-ts`):

| **gene_symbol** | 
| :-------------: | 
| LPL             |   
| FGD5            |
| INPP5B          |

## Output example:

| **gene_symbol** | **Probability** | 
| :-------------: | :--------------:| 
| LPL             |    90.27        |   
| FGD5            |    87.91        |
| INPP5B          |    86.46        |
|   ...           |     ...         |

## License
[MIT](https://choosealicense.com/licenses/mit/)
