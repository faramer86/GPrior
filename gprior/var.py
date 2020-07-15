from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
import pandas as pd
import numpy as np

GTEX_DB = pd.read_csv(
    'databases/gtex_db_v8.tsv', sep='\t')

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
                 'max_features': [None, 'sqrt', 'log2'],
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

NEED = {'DHS',
 'Fantom5',
 'GERP',
 'GTEx',
 'GTEx_Adipose_Subcutaneous',
 'GTEx_Adipose_Visceral_Omentum',
 'GTEx_Adrenal_Gland',
 'GTEx_Artery_Aorta',
 'GTEx_Artery_Coronary',
 'GTEx_Artery_Tibial',
 'GTEx_Brain_Anterior_cingulate_cortex_BA24',
 'GTEx_Brain_Caudate_basal_ganglia',
 'GTEx_Brain_Cerebellar_Hemisphere',
 'GTEx_Brain_Cerebellum',
 'GTEx_Brain_Cortex',
 'GTEx_Brain_Frontal_Cortex_BA9',
 'GTEx_Brain_Hippocampus',
 'GTEx_Brain_Hypothalamus',
 'GTEx_Brain_Nucleus_accumbens_basal_ganglia',
 'GTEx_Brain_Putamen_basal_ganglia',
 'GTEx_Breast_Mammary_Tissue',
 'GTEx_Cells_EBV-transformed_lymphocytes',
 'GTEx_Cells_Transformed_fibroblasts',
 'GTEx_Colon_Sigmoid',
 'GTEx_Colon_Transverse',
 'GTEx_Esophagus_Gastroesophageal_Junction',
 'GTEx_Esophagus_Mucosa',
 'GTEx_Esophagus_Muscularis',
 'GTEx_Heart_Atrial_Appendage',
 'GTEx_Heart_Left_Ventricle',
 'GTEx_Liver',
 'GTEx_Lung',
 'GTEx_Muscle_Skeletal',
 'GTEx_Nerve_Tibial',
 'GTEx_Ovary',
 'GTEx_Pancreas',
 'GTEx_Pituitary',
 'GTEx_Prostate',
 'GTEx_Skin_Not_Sun_Exposed_Suprapubic',
 'GTEx_Skin_Sun_Exposed_Lower_leg',
 'GTEx_Small_Intestine_Terminal_Ileum',
 'GTEx_Spleen',
 'GTEx_Stomach',
 'GTEx_Testis',
 'GTEx_Thyroid',
 'GTEx_Uterus',
 'GTEx_Vagina',
 'GTEx_Whole_Blood',
 'Nearest',
 'PCHiC',
 'Regulome',
 'VEP',
 'VEP_reg',
 'gene_symbol',
 'rank',
 'score',
 'vep_mean',
 'vep_sum'}

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
 'Cells - Cultured fibroblasts': None,
 'Cells - EBV-transformed lymphocytes': None,
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
 'Kidney - Medulla': None,
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
