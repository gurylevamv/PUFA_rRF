import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math

from statannot import add_stat_annotation
from random import sample
from statistics import median
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import f_oneway

from scipy.stats import f_oneway
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multicomp import MultiComparison
from statsmodels.graphics.gofplots import qqplot

import warnings
from IPython.display import display, Math, Latex, Markdown


# # Tumor vs Normal
genes_2cls = pd.read_csv('importance_tumadj_33genes.csv')

test = pd.read_csv('../TEST_adj_tum_TCGA_full.csv', sep='\t')
test = test.set_index('Gene')

ann = ["Tumor" if x[0:3] == "tum" else "Normal" for x in tmp.columns]
pd.Series(ann).value_counts()

test = np.log2(test.T+1)
sns.displot(test.mean())

ind_tumor = [True if x == 'Tumor' else False for x in ann]
ind_adj = [True if x == 'Normal' else False for x in ann]


df_info = pd.DataFrame(columns=['Gene', 'log2FC', 'P-value'])

for gene in genes_2cls['Gene']:
    fold_change = math.log2(test.loc[ind_adj, gene].mean()/test.loc[ind_tumor, gene].mean())
    stat, p = mannwhitneyu(test.loc[ind_adj, gene], test.loc[ind_tumor, gene], alternative='two-sided')
    
    df_info = df_info.append({'Gene': gene, 
                              'log2FC': fold_change, 
                              'P-value': round(p, 3)}, 
                             ignore_index=True)


fdr = fdrcorrection(df_info['P-value'], alpha=0.05)
df_info['P-value_adj'] = fdr[1]


df_info.set_index('Gene',  inplace=True)

for gene in genes_2cls['Gene']:
    ax = sns.violinplot(ann, test[gene],
                     order=['Tumor', 'Normal'], 
                        palette={"Tumor":'#EC3F15', "Normal":'#21DF52'},
                     )

    add_stat_annotation(x=ann, y=test[gene], 
                        order=['Tumor', 'Normal'],
                         box_pairs=[('Tumor', 'Normal')],
                         ax=ax,
                         perform_stat_test=False, pvalues=[df_info.loc[gene, 'P-value_adj']],
                         text_format='star', loc='outside', verbose=0
                       )
    
    
    plt.ylabel(gene, fontsize=24)
    plt.xticks(fontsize=20)
    plt.savefig(fname=f"Adj_Tum/Violin_plot/{gene}.png", format='png')
    plt.close()


# # Molecular Subtypes

genes_4cls = pd.read_csv('../data_for_python/imp_genes_4cls.csv')
genes_4cls.columns = ['Genes', 'Basal', 'Her2', 'LumA', 'LumB', 'MeanDecreaseAccuracy','MeanDecreaseGini']

test_4cls = pd.read_csv('TEST_GSE96058_46_genes.csv',index_col = 0, sep = ',')

lumA, lumB, HER2, Basal = [], [], [], []
for x in test_4cls.index:
    if 'Luminal A' in x:
        lumA.append(x)
    elif 'Luminal B' in x:
        lumB.append(x)
    elif 'Her2' in x:
        HER2.append(x)
    elif 'TNBC' in x:
        Basal.append(x)

subtype = {}
for x in test_4cls.index:
    if 'Luminal A' in x:
        subtype[x] = 'LumA'
    elif 'Luminal B' in x:
        subtype[x] = 'LumB'
    elif 'Her2' in x:
        subtype[x] = 'HER2'
    elif 'TNBC' in x:
        subtype[x] = 'Basal'
        
subtype = pd.Series(subtype)


LumA_ind = subtype[subtype == 'LumA'].index
LumB_ind = subtype[subtype == 'LumB'].index
Her2_ind = subtype[subtype == 'HER2'].index
Basal_ind = subtype[subtype == 'Basal'].index

df_info = pd.DataFrame(columns=['Gene', 'LumA', 'LumB', 'Her2', 'Basal', 'P-value'])

p_vs = []
for gene in genes_4cls.Genes:
    
    means = [test_4cls.loc[LumA_ind, gene].mean(), test_4cls.loc[LumB_ind, gene].mean(), 
             test_4cls.loc[Her2_ind, gene].mean(), test_4cls.loc[Basal_ind, gene].mean()]
    
    stats, p_v = f_oneway(test_4cls.loc[LumA_ind, gene], test_4cls.loc[LumB_ind, gene],
                          test_4cls.loc[Her2_ind, gene], test_4cls.loc[Basal_ind, gene])
    p_vs.append(p_v)
    
    df_info = df_info.append({'Gene': gene, 
                              'LumA': means[0], 
                              'LumB': means[1],
                              'Her2': means[2],
                              'Basal': means[3],
                              'P-value': round(p_v,3)}, 
                             ignore_index=True)


LumA = {}
LumB = {}
Her2 = {}
Basal = {}

for gene in genes_4cls.Genes:
    
    means = [test_4cls.loc[LumA_ind, gene].mean(), test_4cls.loc[LumB_ind, gene].mean(), 
               test_4cls.loc[Her2_ind, gene].mean(), test_4cls.loc[Basal_ind, gene].mean()]
    
    max_mean = max(means)
    
    if means.index(max_mean) == 0:
        pv1 = mannwhitneyu(test_4cls.loc[LumA_ind, gene], test_4cls.loc[LumB_ind, gene], alternative='greater').pvalue
        pv2 = mannwhitneyu(test_4cls.loc[LumA_ind, gene], test_4cls.loc[Her2_ind, gene], alternative='greater').pvalue
        pv3 = mannwhitneyu(test_4cls.loc[LumA_ind, gene], test_4cls.loc[Basal_ind, gene], alternative='greater').pvalue
        
        LumA[gene] = [pv1, pv2, pv3]
            
    elif means.index(max_mean) == 1:
        pv1 = mannwhitneyu(test_4cls.loc[LumB_ind, gene], test_4cls.loc[LumA_ind, gene], alternative='greater').pvalue
        pv2= mannwhitneyu(test_4cls.loc[LumB_ind, gene], test_4cls.loc[Her2_ind, gene], alternative='greater').pvalue
        pv3= mannwhitneyu(test_4cls.loc[LumB_ind, gene], test_4cls.loc[Basal_ind, gene], alternative='greater').pvalue
        
        LumB[gene] = [pv1, pv2, pv3]

    elif means.index(max_mean) == 2:
        pv1 = mannwhitneyu(test_4cls.loc[Her2_ind, gene], test_4cls.loc[LumA_ind, gene], alternative='greater').pvalue
        pv2= mannwhitneyu(test_4cls.loc[Her2_ind, gene], test_4cls.loc[LumB_ind, gene], alternative='greater').pvalue
        pv3= mannwhitneyu(test_4cls.loc[Her2_ind, gene], test_4cls.loc[Basal_ind, gene], alternative='greater').pvalue
        
        Her2[gene] = [pv1, pv2, pv3]
            
    elif means.index(max_mean) == 3:
        pv1 = mannwhitneyu(test_4cls.loc[Basal_ind, gene], test_4cls.loc[LumA_ind, gene], alternative='greater').pvalue
        pv2= mannwhitneyu(test_4cls.loc[Basal_ind, gene], test_4cls.loc[LumB_ind, gene], alternative='greater').pvalue
        pv3= mannwhitneyu(test_4cls.loc[Basal_ind, gene], test_4cls.loc[Her2_ind, gene], alternative='greater').pvalue
        
        Basal[gene] = [pv1, pv2, pv3]


LumA = pd.DataFrame(LumA, index=['p-value1', 'p-value2', 'p-value3']).T
LumB = pd.DataFrame(LumB, index=['p-value1', 'p-value2', 'p-value3']).T
Her2 = pd.DataFrame(Her2, index=['p-value1', 'p-value2', 'p-value3']).T
Basal = pd.DataFrame(Basal, index=['p-value1', 'p-value2', 'p-value3']).T


for df in [LumA, LumB, Her2, Basal]:
    df['p-value1_adj'] = fdrcorrection(df['p-value1'], alpha=0.05)[1]
    df['p-value2_adj'] = fdrcorrection(df['p-value2'], alpha=0.05)[1]
    df['p-value3_adj'] = fdrcorrection(df['p-value3'], alpha=0.05)[1]


LumA = LumA[(LumA['p-value1_adj'] < 0.05) & (LumA['p-value2_adj'] < 0.05) & (LumA['p-value3_adj'] < 0.05)]
LumB = LumB[(LumB['p-value1_adj'] < 0.05) & (LumB['p-value2_adj'] < 0.05) & (LumB['p-value3_adj'] < 0.05)]
Her2 = Her2[(Her2['p-value1_adj'] < 0.05) & (Her2['p-value2_adj'] < 0.05) & (Her2['p-value3_adj'] < 0.05)]
Basal = Basal[(Basal['p-value1_adj'] < 0.05) & (Basal['p-value2_adj'] < 0.05) & (Basal['p-value3_adj'] < 0.05)]


# Create an array with the colors you want to use
colors = ['#55ED81', '#F2FA6C','#F385A8']
colors2= ['#58D68D','#3498DB', '#A569BD', '#E74C3C'] #, '#1ABC9C'
# Set your custom color palette
customPalette = sns.set_palette(sns.color_palette(colors))
customPalette2 = sns.set_palette(sns.color_palette(colors2))


for gene in Basal.index:
    pv1 = mannwhitneyu(test_4cls.loc[Basal_ind, gene], test_4cls.loc[LumA_ind, gene], alternative='greater').pvalue
    pv2 = mannwhitneyu(test_4cls.loc[Basal_ind, gene], test_4cls.loc[LumB_ind, gene], alternative='greater').pvalue
    pv3 = mannwhitneyu(test_4cls.loc[Basal_ind, gene], test_4cls.loc[Her2_ind, gene], alternative='greater').pvalue

    ax = sns.violinplot(subtype, test_4cls[gene], 
                         order=['LumA', 'LumB', 'HER2', 'Basal'], s=1,  xlabel='', stars=False, p_fontsize=8,
                         violin=True,
                         palette=customPalette2, figsize=(7,9), title=f"{gene}", )
    
    add_stat_annotation(x=subtype, y=test_4cls[gene], 
                        order=['LumA', 'LumB', 'HER2', 'Basal'],
                        box_pairs=[( 'Basal', 'LumA'), ( 'Basal', 'LumB'), ('Basal', 'HER2')], 
                        perform_stat_test=False,
                        
                        
                        ax=ax, line_height=0.00, pvalues=[pv1, pv2, pv3],
                        text_format='star', loc='outside', verbose=0
                           )

    plt.xticks(size=10)
    plt.xlabel('Molecular subtype', size=12, labelpad=5)
    plt.tight_layout()
    plt.savefig(fname=f"Violin_subtypes/Basal/{gene}.png", format='png', dpi=120)
    plt.close()