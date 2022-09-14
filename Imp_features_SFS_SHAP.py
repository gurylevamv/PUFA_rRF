import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import math

from statistics import median
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score
from collections import Counter

import shap
import random

random.seed(1)

get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'png'")
plt.rcParams['pdf.fonttype'] = 'truetype'
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.dpi'] = 100
plt.rcParams["figure.figsize"] = (10,7)
sns.set_style("whitegrid", {'axes.grid' : False})

param_grid = { 
    'n_estimators': [450],
    'max_features': ['auto', 'sqrt', 'log2'],
    'max_depth' : [3,5,7,10,15,20,25,50],
    'criterion' :['gini', 'entropy']
}

metabol = pd.read_csv('Pufa_cascade.tsv', sep='\t')

def select_rank(train, test, genes):
    train_sel = train.copy()
    test_sel = test.copy()
    
    train_sel = train_sel.rank()
    test_sel = test_sel.rank()

    train_sel = train_sel.loc[genes,:].T
    test_sel = test_sel.loc[genes,:].T
    
    return train_sel, test_sel

# ## Tumor VS Adjacent

train_2cls = pd.read_csv('../data_for_python/TRAIN_adj_tum_GSE65216_GSE29044_GSE10780.csv' , index_col = 0, sep = '\t')
test_2cls = pd.read_csv('../data_for_python/TEST_adj_tum_TCGA_full.csv',index_col = 0, sep = '\t')

y_train = ["Tumor" if x[0:3]=='tum' else "Adjacent" for x in train_2cls.columns]
y_test = ["Tumor" if x[0:3]=='tum' else "Adjacent" for x in test_2cls.columns]

genes_2cls=pd.read_csv('importance_tumadj_33genes.csv')
genes_2cls.columns = ['Genes', 'Adjacent', 'Tumor', 'MeanDecreaseAccuracy','MeanDecreaseGini']

common = set(train_2cls.index).intersection(test_2cls.index)
train_2cls = train_2cls.loc[common,]
test_2cls = test_2cls.loc[common,]

train_x, test_x = select_rank(train_2cls, test_2cls, genes_2cls.Genes)

# clrf = RandomForestClassifier()
# CV_rf = GridSearchCV(estimator=clrf, param_grid=param_grid, cv=7)
# CV_rf.fit(train_x, y_train)
# best = CV_rf.best_params_

best = {'criterion': 'gini',
 'max_depth': 15,
 'max_features': 'sqrt',
 'n_estimators': 450}


# ## Feature Selection

from mlxtend.feature_selection import SequentialFeatureSelector as SFS
from mlxtend.plotting import plot_sequential_feature_selection as plot_sfs
import matplotlib.pyplot as plt

from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import make_scorer


clrf = RandomForestClassifier(random_state = 42, 
                              max_features = best['max_features'], 
                              n_estimators = 450,
                              max_depth = best['max_depth'], 
                              criterion = best['criterion'])
#clrf.fit(train_x, y_train)

sffs = SFS(clrf, 
           k_features=len(genes_2cls.Genes), 
           forward=True, 
           floating=True, 
           scoring='roc_auc',
           cv=3,
           n_jobs=-1)


sffs = sffs.fit(train_x, y_train)


print('\nSequential Forward Floating Selection (k=3):')
print(sffs.k_feature_idx_)
print('CV Score:')
print(sffs.k_score_)

plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)


fig1 = plot_sfs(sffs.get_metric_dict(), kind='std_dev')

plt.ylim([0.8, 1])
plt.title('Sequential Forward Selection: Tumor vs Adjacent model', fontdict={'fontsize': 18})
plt.ylabel('Performance, ROC-AUC')
plt.grid(alpha=0.3)
plt.show()

# # Molecular Subtypes

train_4cls = pd.read_csv('../data_for_python/TRAIN_GSE25066_GSE81538_GSE31448_GTEx+TCGAnorm.csv' , index_col = 0, sep = ',')
test_4cls = pd.read_csv('../data_for_python/TEST_GSE21653_GSE96058_GTEx+TCGAnorm.csv',index_col = 0, sep = ',')

train_y = []
for x in train_4cls.columns:
    if "Luminal A" in x:
        train_y=train_y+['LumA']
    elif "Luminal B" in x:
        train_y=train_y+['LumB']
    elif "Her2" in x:
        train_y=train_y+['HER2']
    else:
        train_y=train_y+['Basal']        

test_y = []
for x in test_4cls.columns:
    if "Luminal A" in x:
        test_y=test_y+['LumA']
    elif "Luminal B" in x:
        test_y=test_y+['LumB']
    elif "Her2" in x:
        test_y=test_y+['HER2']
    else:
        test_y=test_y+['Basal'] 


genes_4cls=pd.read_csv('imp_genes_4cls.csv')
genes_4cls.columns = ['Genes', 'Basal', 'Her2', 'LumA', 'LumB', 'MeanDecreaseAccuracy','MeanDecreaseGini']

common = set(train_4cls.index).intersection(test_4cls.index)
train_4cls = train_4cls.loc[common,]
test_4cls = test_4cls.loc[common,]


train_x, test_x = select_rank(train_4cls, test_4cls, genes_4cls.Genes)


# clrf = RandomForestClassifier()
# CV_rf = GridSearchCV(estimator=clrf, param_grid=param_grid, cv=7)
# CV_rf.fit(train_x, train_y)
# best = CV_rf.best_params_ 

best = {'criterion': 'gini',
 'max_depth': 50,
 'max_features': 'sqrt',
 'n_estimators': 450}


# ## SHAP values molecular Subtypes

clrf = RandomForestClassifier(random_state = 1, 
                              max_features = best['max_features'], 
                              n_estimators = 450,
                              max_depth = best['max_depth'], 
                              criterion = best['criterion'])
clrf.fit(train_x, train_y)


explainer = shap.TreeExplainer(clrf)
shap_values = explainer.shap_values(train_x, approximate=False, check_additivity=False)

shap.summary_plot(shap_values,  train_x, class_names=clrf.classes_,plot_type="bar")


#Basal
shap.summary_plot(shap_values[0], train_x, max_display=5, plot_type="dot")
#HER2
shap.summary_plot(shap_values[1], train_x, max_display=5, plot_type="dot")
#LumA
shap.summary_plot(shap_values[2], train_x, max_display=5, plot_type="dot")
#LumB
shap.summary_plot(shap_values[3], train_x, max_display=5, plot_type="dot")

