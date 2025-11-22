from Bio.Seq import Seq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
#Creating Python Function To Convert DNA into Protein
def convert_dna_to_protein(dna_sequence):

    #Ensuring the DNA sequence is valid
    for i in dna_sequence:
        if i != "A" and i != "T" and i != "G" and i != "C":
            raise ValueError("Invalid DNA sequence")
            quit()
        else:
            continue

    #Converting DNA sequence to Protein sequence
    protein_sequence = Seq(dna_sequence).translate()
    return protein_sequence

#dna_sequence = "ATGCTGCCCBTAAAA"
#protein_sequence = convert_dna_to_protein(dna_sequence)
#print(f"Protein Sequence: {protein_sequence}")

#Crating Function to calculate hamming distance between slack username and twitter handle
def calculate_hamming_distance(slack_username, twitter_handle):
    hamming_distance = 0
    if len(slack_username) != len(twitter_handle):
        raise ValueError("Slack username and Twitter handle must have the same length")
    else:
        for i in range(len(slack_username)):
            if slack_username[i] != twitter_handle[i]:
                hamming_distance += 1
        return hamming_distance
#Example usage
#slack_username = "@pavanpunnamraju"
#twitter_handle = "pavan_punnamraju"
#hamming_distance = calculate_hamming_distance(slack_username, twitter_handle)
#print(f"Hamming Distance: {hamming_distance}")

#Initializing Plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

##Part A - Gene Expression Analysis

#Loading Gene Expression Data
df = pd.read_csv("hbr_uhr_top_deg_normalized_counts.csv")
#print(gene_expression_data.head())

#Setting the index of the dataframe to the first column
gene_expression_data = df.set_index("Unnamed: 0")
#plotting clustured heatmap of hbr vs uhr gene expression data
sns.clustermap(gene_expression_data, linewidths=0.5, linecolor="black", cmap="Blues")

#--------------------------------------------------------------------------------------------------------
##Part A2
#--------------------------------------------------------------------------------------------------------

#Loading DEG data
df_two = pd.read_csv("hbr_uhr_deg_chr22_with_significance.csv")
palette = {"down":"orange", "up":"green", "ns":"grey"}
sns.scatterplot(x="log2FoldChange", y="-log10PAdj", hue = "significance", palette = palette, data=df_two, ax=axes[0,1])

#Dashed lines
axes[0,1].axhline(y=0, color="grey", linestyle="--")
axes[0,1].axvline(x=1, color="grey", linestyle="--")
axes[0,1].axvline(x=-1, color = 'grey', linestyle = "--")

#General formatting
axes[0,1].set_xlabel("log2FoldChange")
axes[0,1].set_ylabel("-log10PAdj")
axes[0,1].legend(title="Significance")

#-------------------------------------------------------------------------------------------------------
##Part B1
#-------------------------------------------------------------------------------------------------------
#Loading dataset
brst_cancer_df = pd.read_csv('data-3.csv')
#print(brst_cancer_df.columns)

#Setting palette
palette = {"M":'blue', "B":'orange'}

#Plotting data
sns.scatterplot(x="radius_mean", y="texture_mean", hue = "diagnosis", palette = palette, ax = axes[0,2],data = brst_cancer_df)

#-------------------------------------------------------------------------------------------------------
##Part B2
#-------------------------------------------------------------------------------------------------------
brst_cancer_matrix = brst_cancer_df[['radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean']]
correlation_matrix = brst_cancer_matrix.corr()
sns.heatmap(data = correlation_matrix, cmap = "Blues",linecolor='black', linewidths=0.5, annot = True, ax = axes[1,0])


#-------------------------------------------------------------------------------------------------------
##Part B3
#-------------------------------------------------------------------------------------------------------

#Setting palette
palette = {"M":'blue', "B":'orange'}

#Plotting data
sns.scatterplot(x="smoothness_mean", y="compactness_mean", hue = "diagnosis", palette = palette, ax = axes[1,1], data = brst_cancer_df)


#-------------------------------------------------------------------------------------------------------
##Part B4
#-------------------------------------------------------------------------------------------------------

#Setting palette
palette = {"M":'blue', "B":'orange'}

#Converting area mean into gaussian distribution
density = gaussian_kde(brst_cancer_df['area_mean'])

#Plotting data
sns.kdeplot(x='area_mean', hue = 'diagnosis', palette = palette, ax = axes[1,2], fill = True, data = brst_cancer_df)

fig.tight_layout()
plt.show()
