from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
import pandas as pd
import numpy as np

GTEX_DB = pd.read_csv(
    'databases/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.tsv', sep='\t')

GTEX_SIMILARITY_DB = pd.read_csv(
    'databases/UCSC_gtex_db.tsv', sep='\t', index_col=0)

BLASTP_SIMILARITY_DB = pd.read_csv(
    'databases/UCSC_blastp_db.tsv', sep='\t', index_col=0)

ATLAS_SIMILARITY_DB = pd.read_csv(
    'databases/UCSC_atlas_db.tsv', sep='\t', index_col=0)

GENE_INTERACTIONS_DB = pd.read_csv(
    'databases/UCSC_gene_interactions.tsv', sep='\t', index_col=0)

REACTOME_DB = pd.read_csv(
    'databases/FIsInGene_122718_with_annotations.txt', sep='\t')

PARAM_DIST_RF = {'max_depth': [5, 10, 30, 50],
                 'max_features': ['auto', 'sqrt', 'log2'],
                 'min_samples_split': [2, 5],
                 'min_samples_leaf': [2, 3]}


PARAM_DIST_SVC = {'gamma': [1, 0.1, 0.01, 0.001, 0.0001],
                  'C': [0.001, 0.01, 0.1, 1]}

PARAM_DIST_LR = {"C": np.logspace(-5, 3, 40),
                 "penalty": ["l1", "l2"]}

PARAM_DIST_DT = {'max_depth': [5, 10, 15, 50, 100],
                 'min_samples_leaf': [1, 2, 3, 4],
                 'min_samples_split': [2, 3, 5, 10]}

PARAM_DIST_ADA = {"learning_rate": np.logspace(-6, 0, 15),
                  "algorithm": ['SAMME', 'SAMME.R']}

WEIGHTS = {1:9, 0:1}
 
THRESHOLD_RANGE = [0.01, 0.1] + [i for i in range(1, 98)]

PARAMS = {
    
    'Support Vector Machine': PARAM_DIST_SVC,
    'ADABoosting': PARAM_DIST_ADA,
    'Random Forest': PARAM_DIST_RF,
    'Logistic regression': PARAM_DIST_LR,
    'Decision Tree' : PARAM_DIST_DT 
}

ADA_BASE = DecisionTreeClassifier(criterion='gini', class_weight=WEIGHTS)

MODELS = {
    'Logistic regression': [LogisticRegression(solver='liblinear', class_weight=WEIGHTS), PARAM_DIST_LR],    
    'Support Vector Machine': [SVC(probability=True, C=0.01, gamma=0.0001, class_weight=WEIGHTS), PARAM_DIST_SVC],
    'ADABoosting': [AdaBoostClassifier(base_estimator = ADA_BASE, learning_rate=0.01), PARAM_DIST_ADA],  
    'Random Forest': [RandomForestClassifier(bootstrap=True, n_estimators=10, class_weight=WEIGHTS, oob_score=False), PARAM_DIST_RF],
    'Decision Tree': [DecisionTreeClassifier(class_weight=WEIGHTS), PARAM_DIST_DT]
}

DO_NOT_NEED = ['ld_snp_rsID',
               'chrom',
               'pos',
               'GRCh38_chrom',
               'gwas_size',
               'GRCh38_pos',
               'GRCh38_gene_chrom',
               'gene_chrom',
               'gene_tss',
               'GRCh38_gene_pos',
               'disease_name',
               'disease_efo_id',
               'cluster_id',
               'gwas_source',
               'gwas_snp',
               'gwas_pvalue_description',
               'gwas_odds_ratio',
               'gwas_odds_ratio_ci_start',
               'gwas_odds_ratio_ci_end',
               'gwas_size',
               'gwas_pmid',
               'gwas_study',
               'gwas_reported_trait',
               'vep_terms',
               'gnomad',
               'gnomad_sas',
               'gnomad_oth',
               'gnomad_asj',
               'gnomad_nfe',
               'gnomad_afr',
               'gnomad_amr',
               'gnomad_fin',
               'gnomad_eas',
               'afr',
               'amr',
               'eas',
               'eur',
               'sas',
               'gwas_beta',
               'gwas_pvalue',
               'ls_snp_is_gwas_snp',
               'r2']

AGG = {
 'DHS':'max',
 'Fantom5':'max',
 'GERP':'max',
 'GTEx':'max',
 'GTEx_Adipose_Subcutaneous':'max',
 'GTEx_Adipose_Visceral_Omentum':'max',
 'GTEx_Adrenal_Gland':'max',
 'GTEx_Artery_Aorta':'max',
 'GTEx_Artery_Coronary':'max',
 'GTEx_Artery_Tibial':'max',
 'GTEx_Brain_Anterior_cingulate_cortex_BA24':'max',
 'GTEx_Brain_Caudate_basal_ganglia':'max',
 'GTEx_Brain_Cerebellar_Hemisphere':'max',
 'GTEx_Brain_Cerebellum':'max',
 'GTEx_Brain_Cortex':'max',
 'GTEx_Brain_Frontal_Cortex_BA9':'max',
 'GTEx_Brain_Hippocampus':'max',
 'GTEx_Brain_Hypothalamus':'max',
 'GTEx_Brain_Nucleus_accumbens_basal_ganglia':'max',
 'GTEx_Brain_Putamen_basal_ganglia':'max',
 'GTEx_Breast_Mammary_Tissue':'max',
 'GTEx_Cells_EBV-transformed_lymphocytes':'max',
 'GTEx_Cells_Transformed_fibroblasts':'max',
 'GTEx_Colon_Sigmoid':'max',
 'GTEx_Colon_Transverse':'max',
 'GTEx_Esophagus_Gastroesophageal_Junction':'max',
 'GTEx_Esophagus_Mucosa':'max',
 'GTEx_Esophagus_Muscularis':'max',
 'GTEx_Heart_Atrial_Appendage':'max',
 'GTEx_Heart_Left_Ventricle':'max',
 'GTEx_Liver':'max',
 'GTEx_Lung':'max',
 'GTEx_Muscle_Skeletal':'max',
 'GTEx_Nerve_Tibial':'max',
 'GTEx_Ovary':'max',
 'GTEx_Pancreas':'max',
 'GTEx_Pituitary':'max',
 'GTEx_Prostate':'max',
 'GTEx_Skin_Not_Sun_Exposed_Suprapubic':'max',
 'GTEx_Skin_Sun_Exposed_Lower_leg':'max',
 'GTEx_Small_Intestine_Terminal_Ileum':'max',
 'GTEx_Spleen':'max',
 'GTEx_Stomach':'max',
 'GTEx_Testis':'max',
 'GTEx_Thyroid':'max',
 'GTEx_Uterus':'max',
 'GTEx_Vagina':'max',
 'GTEx_Whole_Blood':'max',
 'Nearest':'max',
 'PCHiC':'max',
 'Regulome':'max',
 'VEP':'max',
 'VEP_reg':'max',
 'rank':'min',
 'score':'max',
 'vep_mean':'max',
 'vep_sum':'max',
}

AGG_mean = {
 'DHS':'mean',
 'Fantom5':'mean',
 'GERP':'mean',
 'GTEx':'mean',
 'Nearest':'mean',
 'PCHiC':'mean',
 'Regulome':'mean',
 'VEP':'mean',
 'VEP_reg':'mean',
 'rank':'mean',
 'score':'mean'} 

AGG_mean_names = ['DHS_mean',
 'Fantom5_mean',
 'GERP_mean',
 'GTEx_mean',
 'Nearest_mean',
 'PCHiC_mean',
 'Regulome_mean',
 'VEP_mean',
 'VEP_reg_mean',
 'rank_mean',
 'score_mean']

GTEX_COLUMNS = {'Adipose - Subcutaneous': None,
                'Adipose - Visceral (Omentum)': None,
                'Adrenal Gland': None,
                'Artery - Aorta': None,
                'Artery - Coronary': None,
                'Artery - Tibial': None,
                'Bladder': None,
                'Brain - Amygdala': None,
                'Brain - Anterior cingulate cortex (BA24)': None,
                'Brain - Caudate (basal ganglia)': None,
                'Brain - Cerebellar Hemisphere': None,
                'Brain - Cerebellum': None,
                'Brain - Cortex': None,
                'Brain - Frontal Cortex (BA9)': None,
                'Brain - Hippocampus': None,
                'Brain - Hypothalamus': None,
                'Brain - Nucleus accumbens (basal ganglia)': None,
                'Brain - Putamen (basal ganglia)': None,
                'Brain - Spinal cord (cervical c-1)': None,
                'Brain - Substantia nigra': None,
                'Breast - Mammary Tissue': None,
                'Cells - EBV-transformed lymphocytes': None,
                'Cells - Transformed fibroblasts': None,
                'Cervix - Ectocervix': None,
                'Cervix - Endocervix': None,
                'Colon - Sigmoid': None,
                'Colon - Transverse': None,
                'Esophagus - Gastroesophageal Junction': None,
                'Esophagus - Mucosa': None,
                'Esophagus - Muscularis': None,
                'Fallopian Tube': None,
                'Heart - Atrial Appendage': None,
                'Heart - Left Ventricle': None,
                'Kidney - Cortex': None,
                'Liver': None,
                'Lung': None,
                'Minor Salivary Gland': None,
                'Muscle - Skeletal': None,
                'Nerve - Tibial': None,
                'Ovary': None,
                'Pancreas': None,
                'Pituitary': None,
                'Prostate': None,
                'Skin - Not Sun Exposed (Suprapubic)': None,
                'Skin - Sun Exposed (Lower leg)': None,
                'Small Intestine - Terminal Ileum': None,
                'Spleen': None,
                'Stomach': None,
                'Testis': None,
                'Thyroid': None,
                'Uterus': None,
                'Vagina': None,
                'Whole Blood': None}