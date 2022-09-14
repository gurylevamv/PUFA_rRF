#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt


# In[ ]:


from gseapy.plot import barplot, dotplot


# In[ ]:


exlude_4cls = ["CYP2S1","CYP4A22" ,
"CYP4F22",
"CYP4V2",
"CYP4X1",
'CYP4Z1',
'CYP3A7',
'CYP21A2',
'CYP26C1',
'CYP27C1' ,
'FAR1',
'ACOT1',
'ACOT2',
'ACOT4',
'ACOT6',
'ACOT12',
'FADS1',
'FADS6',
'ACAD9',
'ACAD11',
'PRXL2B',
'PTGR1',
'LTA4H',
'EPHX4',
'FAAH2',
'PLA2G2C',
'PLA2G4B',
'PLA2G4D',
'PLA2G4E',
'PLA2G4F',
'PLA2G10',
'PLA2G12B',
'PLB1',
'PLCZ1',
'PLCD3',
'PLCD4',
'PLD4',
'PLD5',
'PLD6',
'GPR55',
'LGR6',
'FFAR1',
'FFAR4',
'OXER1',
'OXGR1',
'THEM4',
'THEM5']

PUFA = pd.read_csv('PUFA_gene_list.tsv', sep='\t', index_col=0)
backgoround_4cls= PUFA[~PUFA.Gene.isin(exlude_4cls)].Gene


# In[2]:


LumA = ['ELOVL5',
 'ACAA1',
 'PLD2',
 'ACAD8',
 'PLCL1',
 'CYP4F11',
 'PTGER3',
 'CYP4F8',
 'ELOVL2',
 'EPHX2',
 'LPCAT3',
 'LTC4S']

LumB = ['PTGES3', 'ADIPOR1', 'MBOAT7', 'PLAA', 'ACOT7', 'ACOT8', 'CYP2B6', 'FAAH']

Her2 = ['FASN', 'FABP6', 'PLCH1', 'PLCB4', 'ALOX15B', 'FADS2']

Basal = ['AKR1B1',
 'CYP39A1',
 'PLD1',
 'PLA2G4A',
 'FPR2',
 'PLCG2',
 'CYP1B1',
 'CYP7B1',
 'FABP5',
 'PLA2G7',
 'CBR1']


# In[ ]:


libraries = ['KEGG_2021_Human', 
                            'BioCarta_2016', 
                            'WikiPathway_2021_Human',
                            'Panther_2016',
                            'GO_Biological_Process_2021', 
                            'GO_Molecular_Function_2021']


# In[ ]:


enr_4cls = gp.enrichr(gene_list=Basal, #select subtype gene list
#                 description='test_name',
                 gene_sets=libraries,
                 background=backgoround_4cls, # or the number of genes, e.g 20000
                 outdir='LumA/',
                 cutoff=0.001, # only used for testing.
                 format='png',
                 verbose=False)


# In[ ]:


#Basal
for geneset in libraries:
    if not enr_4cls.results[enr_4cls.results.Gene_set==geneset].empty:
        ax = dotplot(enr_4cls.results[enr_4cls.results.Gene_set==geneset], title=geneset, cmap='viridis_r')
        #ax.grid(False)
#         plt.savefig(f"test/test2/Basal/Basal_{geneset}.png", format='png')
#         plt.close()
        plt.show()


# In[ ]:


#Her2
for geneset in libraries:
    if not enr_4cls.results[enr_4cls.results.Gene_set==geneset].empty:
        ax = dotplot(enr_4cls.results[enr_4cls.results.Gene_set==geneset], title=geneset, cmap='viridis_r')
        #ax.grid(False)
#         plt.savefig(f"test/test2/Her2/Her2_{geneset}.png", format='png')
#         plt.close()
        plt.show()


# In[ ]:


#LumB
for geneset in libraries:
    if not enr_4cls.results[enr_4cls.results.Gene_set==geneset].empty:
        ax = dotplot(enr_4cls.results[enr_4cls.results.Gene_set==geneset], title=geneset, cmap='viridis_r')
        #ax.grid(False)
#         plt.savefig(f"test/test2/LumB/LumB_{geneset}.png", format='png')
#         plt.close()
        plt.show()


# In[ ]:


#LumA
for geneset in libraries:
    if not enr_4cls.results[enr_4cls.results.Gene_set==geneset].empty:
        ax = dotplot(enr_4cls.results[enr_4cls.results.Gene_set==geneset], title=geneset, cmap='viridis_r')
        #ax.grid(False)
#         plt.savefig(f"test/test2/LumA/LumA_{geneset}.png", format='png')
#         plt.close()
        plt.show()


# # Tumor vs Normal

# In[ ]:


exclude_2cls=[
'CYP4A22',
'CYP3A7',
'CYP21A2',
'CYP26C1',
'CYP27C1',
'ACOT1',
'ACOT2',
'FADS1',
'FADS6',
'ACAD11',
'PRXL2B',
'LTA4H',
'PLA2G2C',
'PLA2G4B',
'PLA2G4E',
'PLA2G10',
'FFAR4'
]

backgoround_2cls= PUFA[~PUFA.Gene.isin(exclude_2cls)].Gene


# In[ ]:


genes_up_in_tumor = pd.read_csv('UP_in_Tumor_vs_Adj.tsv', sep='\t', index_col=0)


# In[ ]:


enr2 = gp.enrichr(gene_list=genes_up_in_tumor.Gene,
                 # or gene_list=glist
                 description='test_name',
                 gene_sets=libraries,
                 background=backgoround_2cls, # or the number of genes, e.g 20000
                 outdir='Enrichr_Tumor/',
                 cutoff=0.01, # only used for testing.
                 format='png',
                 verbose=False)


# In[ ]:


for geneset in libraries:
    if not enr2.results[enr2.results.Gene_set==geneset].empty:
        ax = dotplot(enr2.results[enr2.results.Gene_set==geneset], title=geneset, cmap='viridis_r')
        #ax.grid(False)
#         plt.savefig(f"Adj_Tum_BRCA/Enrichr_Adj/UPreg_Adj_{geneset}.png", format='png')
#         plt.close()
        plt.show()


# In[ ]:


genes_up_in_adj = pd.read_csv('UP_in_Adj_vs_Tumor.tsv', sep='\t', index_col=0)


# In[ ]:


enr2 = gp.enrichr(gene_list=genes_up_in_adj.Gene,
                 # or gene_list=glist
                 description='test_name',
                 gene_sets=libraries,
                 background=backgoround_2cls, # or the number of genes, e.g 20000
                 outdir='Enrichr_Adj/',
                 cutoff=0.01, # only used for testing.
                 format='png',
                 verbose=False)


# In[ ]:


for geneset in libraries:
    if not enr2.results[enr2.results.Gene_set==geneset].empty:
        ax = dotplot(enr2.results[enr2.results.Gene_set==geneset], title=geneset, cmap='viridis_r')
        plt.show()

