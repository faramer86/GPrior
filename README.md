![](logo.png)

# **Gene Prioritization tool (GPrior)**

GPrior is a Gene Prioritization tool that uses PU learning to prioritize candidate genes based on their similarity with a set of known positive examples. As a result it returns a table with probabilities for each gene.

## **Table of Contents**

1. [Dependencies](#dependencies)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Features](#features)
5. [Input example](#input-example)
6. [Output example](#output-example)
7. [References](#reference)

## **Dependencies** 

Make sure that you have `python` version >= 3.6 and latest version of `pip`. You can check it by running:

```bash
 python --version
 pip --version
```

If you do not have Python, please install the latest 3.x version from python.org.

In order to be cross-platform we uses `virtualenv`. If you do not have it you can install it using `pip` (use `sudo` if necessary):

```bash
pip install virtualenv
```

If you have any problems with the installation or usage of `virtualenv` be free to consult the official documentation: https://virtualenv.pypa.io/en/latest/

## **Installation**

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

**Now you can launch `gprior.py`, it should work!**

#### NOTE:

1) Do not forget to exit virtual environment after you are done:

```bash
deactivate
```

2) Do not forget to activate `vprior` everytime you want to launch tool:

```bash
source vprior/bin/activate
```

## **Usage**

GPrior is a commandline tool that provides user interface for tuning and training Positive-Unlabeled classifier. Minimum input: table with positive examples (gene symbols) and dataframe (rows - gene symbols, columns - gene features). For more details/examples see sections below. In order to perform prioritization, user can try to repeat our pipeline from the article (see **Reference** section) or compile his own table with features. 

There are two main scripts: `gprior.py` and `process_postgap.py`. If you want to prioritize genes based on your own table of features - use `gprior.py`. If you want to reproduce pipeline from the article - use `process_postgap.py` for postgap output processing and only then use `gprior.py`. For more details read futher.

### process_postgap.py

This script allows you to summarize postgap output and make gene-level table for subsequent analysis. All the arguments are pretty straigtforward. For more details see **Input examples** section.

```bash
  -h, --help      show this help message and exit
  -i , --input    Path to file/folder with postgap output
  -o , --output   Path to output
```

### gprior.py

`gprior.py` perfor gene prioritization using table of features and sets of causal genes as an input (see **Input examples** section). 

It consists of 5 PU Bagging classifiers. All the predictions are combined using optimal combination approach (see article) or simple mean. The main difference between this two approaches is that if you do not have enough data for AES compiling, so you can't evaluate each PU bagging prediction quality and find the best combination, so then simple mean between all the predictions is calculated (just do not specify `-aes` argument). Otherwise, true gene set (TS) is used for training and algorithm selection set (ASS) for optimal combination finding and quality evaluation.

We add several useful tunable parameters. See `--help` page:

### Arguments:

```bash
  -h, --help            show this help message and exit
  -i , --input          Path to gene-features table
  -ts , --true_set      Path to table with gene symbols of causal genes (True
                        gene set)
  -aes , --algorithm_evaluation_set 
                        Path to file with algorithm evaluation set
  -o , --output         Path to output file
  -n , --n_bootstrap    (default=15) Number of bootstraps for PU Bagging
  -k , --k_clusters     (default=n_features) Number of clusters for Feature
                        Agglomeration
  -s , --s_coef         (default=1) Sampling coefficient
  --add_features        (positional; default=False) Add additional featutes to
                        the provided table
  --tune                (positional; default=False) Tune hyperparameters of
                        the model It significently slow down training process.
                        Do not use it with high n.
  --set_seed            (positional; default=False) Switch it on if you want
                        reproducibility
```
example:

1) launching **without** AES:

```bash
 ./gprior.py \
    -i test_input_copy.tsv \ 
    -ts test_ts.tsv \ 
    -o test_output.tsv 
```

2) launching **with** AES:

```bash
 ./gprior.py \ 
    -i test_input_copy.tsv \ 
    -ts test_ts.tsv \ 
    -aes test_ass.tsv \ 
    -o test_output.tsv \
    -n 50
```

## **Input example:**

There are two main requirement for the input files:

1) **All files have to have .tsv extansion**

2) **All files have to have column with gene symbols, called "gene_symbol"**

* `./gprior.py` Input file (`-i`):

 |**gene_symbol**| **feature 1** |**feature 2**|**...**|**feature n**|
 |:----:| :--------------------: |:--------------------:|---|:--------------------:|
 |GENE SYMBOL 1| ... </br>  | ... </br> |...|... </br> |
 |...|...| ... |...|...|...|
 |GENE SYMBOL n| ... </br> |... </br> |...|... </br>|

* `./process_postgap.py` Input file (`-i`):

For more detailes see: https://github.com/Ensembl/postgap

 |**gene_symbol**| ** SNPs ** |**feature 2**|**...**|**feature n**|
 |:----:| :--------------------: |:--------------------:|---|:--------------------:|
 |GENE SYMBOL 1| SNP1 </br>  | ... </br> |...|... </br> |
 |GENE SYMBOL 1| SNP2 </br>  | ... </br> |...|... </br> |
 |GENE SYMBOL 1| SNP3 </br>  | ... </br> |...|... </br> |
 |...|...| ... |...|...|...|
 |GENE SYMBOL n| ... </br> |... </br> |...|... </br>|

* Algorithm selection set (`-ass`) and True gene set (`-ts`):

| **gene_symbol** | 
| :-------------: | 
| LPL             |   
| FGD5            |
| INPP5B          |

## **Output example:**

| **gene_symbol** | **Probability** | 
| :-------------: | :--------------:| 
| LPL             |    90.27        |   
| FGD5            |    87.91        |
| INPP5B          |    86.46        |
|   ...           |     ...         |

## **License**
[MIT](https://choosealicense.com/licenses/mit/)

## **Reference**

-//-
