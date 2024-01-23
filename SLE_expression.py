import numpy as np
import pandas as pd
import rpy2.robjects as ro
import seaborn as sns
import matplotlib.pyplot as plt
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from sklearn.decomposition import PCA
path_to_limma = "/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library" #importing limma library from R
limma = importr('limma', lib_loc=path_to_limma) #importing the writexl library from R
writexl = importr('writexl', lib_loc="/Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library")

GSE17755 = pd.read_csv('GSE17755_series_matrix.txt', delimiter='\t', skiprows = 34, low_memory=False) #getting the gene expression data

#This matrix has other autoimmune patients in it so needs to be cleaned
colnames = GSE17755.columns 
save_GSE17755 = ['!Sample_title']
control_vs_SLE = {} #dictionary to determine if sample is control or SLE
for col in colnames:
    if col.startswith("SLE"):
        save_GSE17755.append(col)
        control_vs_SLE[GSE17755[col].iloc[0]] = "SLE"
    elif col.startswith("Control healthy individual"):
        save_GSE17755.append(col)
        control_vs_SLE[GSE17755[col].iloc[0]] = "Control"
GSE17755 = GSE17755[save_GSE17755]

# Setting rownames
GSE17755 = GSE17755[GSE17755['!Sample_title'].str.startswith('AGhs')]
GSE17755 = GSE17755.dropna()
GSE17755 = GSE17755.set_index(GSE17755.columns[0]) #has SLE (n=22), Control (n=45)
# Log2 ratios of Cy3 to Cy5 were calculated and normalized by the method of global ratio median normalization

#Getting data for probe to gene for GSE17755
probe_info_starts = 320
probe_info_ends = 30657
GSE17755_probe_to_gene = pd.read_csv('GSE17755_family.csv', delimiter='\t', skiprows = 319, nrows = 30336)

pca_df = GSE17755.T #transposing data so samples are on rows and probes are on columns

#Extracting condition for PCA plot
rownames = pca_df.index.tolist()
sample_type_pca = []
for row in rownames:
    if row.startswith("SLE"):
        sample_type_pca.append("SLE")
    elif row.startswith("Control"):
        sample_type_pca.append("Control")
        
pca = PCA(n_components=2) #PC1 and PC2
pc_scores = pca.fit_transform(pca_df)

pc_df = pd.DataFrame(data={'PC1': pc_scores[:, 0], 'PC2': pc_scores[:, 1], 'Condition': sample_type_pca})

#Variance explained by PCs
var_explained = pca.explained_variance_ratio_

#PCA Plot
sns.set(style='whitegrid')
plt.figure(figsize=(6, 4))
sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pc_df, palette='colorblind', s=100)
#Adding percentage of variance explained to plot
plt.xlabel(f'PC1 ({var_explained[0]*100:.2f}%)')
plt.ylabel(f'PC2 ({var_explained[1]*100:.2f}%)')
plt.title('Variance Explained by 1st and 2nd Components for Normalized Data')
plt.legend()
plt.show()
# Difference between disease and control loose clustering, can say overlap due to housekeeping genes
# Metric to look at how close clusters are (gives two clusters and calculates)

#When |log2(FC)|> 1.5 and p-value < 0.05, the gene is considered as being statistically significant in the 
#combined matrix
GSE17755 #samples on cols and probes as rows (correct format for limma)
sample_names = GSE17755.columns
labels = []
sample_ids = []
#Must make sample_ids numeric 
for i in range(len(sample_names)):
    if sample_names[i].startswith("SLE"):
        label = 1  #SLE
        sample_id = 1+float(sample_names[i].split(" ")[2])/100 #Unique sample ID from SLE
    if sample_names[i].startswith("Control"):
        label = 0  #Control
        sample_id = float(sample_names[i].split(" ")[3])/100  #Unique sample ID from controls
    labels.append(label)
    sample_ids.append(sample_id)

#Making the design matrix for limma
design = pd.DataFrame({'Label': labels, 'Sample_ID': sample_ids})
GSE17755.columns = list(design["Sample_ID"])
GSE17755 = GSE17755.apply(pd.to_numeric, errors='coerce') #ensuring values are input to R as numbers

with (ro.default_converter + pandas2ri.converter).context():
    r_data = ro.conversion.get_conversion().py2rpy(GSE17755) #making GSE17755 data into R formating
    r_design = ro.conversion.get_conversion().py2rpy(design) #making design matrix into R formating
    
fit = limma.lmFit(r_data, r_design) #fitting data to linear model according to design

contrast_matrix = limma.makeContrasts("1 - 0", levels=r_design) #since SLE is 1 and 0 is control
fit2 = limma.lmFit(r_data, r_design)

fit2 = limma.eBayes(fit2) #Bayesian adjustment

r_output = limma.topTreat(fit2, coef=1, number=np.Inf)
writexl.write_xlsx(r_output, "limma_output.xlsx") #writing to an xlsx file so can be re-uploaded to python 

# Since excel loses the rownames need to get them as a list.
row_names = r_output.rownames
row_names_list = list(row_names)

limma_data = pd.read_csv('limma_output.csv') #reading in limma output to pandas dataframe

limma_data.index = row_names_list #setting probe names since they were not included in excel file
#explain logic behind paper deviations

#Setting Up a Column to indicate if the gene is significant or not
limma_data['Change'] = 'Not Sig'
limma_data.loc[(limma_data['logFC'] > 1.5) & (limma_data['adj.P.Val'] < 0.05), 'Change'] = 'Up'
limma_data.loc[(limma_data['logFC'] < -1.5) & (limma_data['adj.P.Val'] < 0.05), 'Change'] = 'Down'
limma_data['-log10(pvalue)'] = -np.log10(limma_data['adj.P.Val'])

#Volcano Plot
sns.set(style='whitegrid')
plt.figure(figsize=(6, 4))
sns.scatterplot(x='logFC', y='-log10(pvalue)', hue='Change', data=limma_data, palette='colorblind', s=5)
plt.xlabel('log2(Fold Change)')
plt.ylabel('-log10(p-value)')
plt.title('Differential Gene Expression')
plt.legend()
plt.show()

#Taking out Only Significant Probes
filtered_limma = limma_data[(limma_data['adj.P.Val'] < 0.05) & (abs(limma_data['logFC']) > 1.5)]
probes_sig = list(filtered_limma.index)
#Getting the metadata for significant probes
GSE17755_probe_to_gene
filtered_probe_to_gene = GSE17755_probe_to_gene[GSE17755_probe_to_gene['ID'].isin(probes_sig)]
#Dataframe of the relevant gene information for the significant genes found
sig_gene_info = filtered_probe_to_gene[["ID", "GeneName", "Entrez Gene ID", "Symbol", 
                                        "GO Molecular Function", "GO Biological Process", 
                                        "GO Cellular Component"]]

result_df = pd.merge(sig_gene_info, filtered_limma, left_on='ID', right_index=True, how='inner')
#getting metadata and significant limma results into same dataframe

#Printing out the top 10 LogFC genes identified through limma
results_sorted = result_df.reindex(result_df['logFC'].abs().sort_values(ascending=False).index)
results_sorted = results_sorted.dropna(subset=['Symbol'])
top_10 = results_sorted.head(10)
display_table = top_10[["Symbol", "GeneName", "logFC"]]
display_table.to_csv("top_10_DE_genes.csv", index = False)

def get_go_terms_counts(col_name):
    #creating a function to find out what pathways the significantly expressed genes belong to 
    split_df = result_df[col_name].str.split(', ', expand=True)
    melted_df = split_df.melt(value_name='value').dropna()['value']
    
    go_counts = {}
    for i in melted_df:
        if i == "nan":
            continue
        if i in go_counts.keys():
            go_counts[i] += 1
        else:
            go_counts[i] = 1
    
    go_counts = dict(sorted(go_counts.items(), key=lambda item: item[1], reverse=True))
    
    return(go_counts)

go_mf_counts = get_go_terms_counts("GO Molecular Function")
go_bp_counts = get_go_terms_counts("GO Biological Process")
go_cell_components_counts = get_go_terms_counts("GO Cellular Component")