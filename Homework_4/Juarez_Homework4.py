import pandas as pd

'''
INSERT CODE HERE
Homework 3.A
(i) Merge the 5 normal CSV files together and the 5 tumor CSV files, 
result should 2 separate Dataframes, one with Normal variants and another with Tumor variants.

SUGGESTION: Prior to coding, create 2 empty folders, Normal_CSV and Tumor_CSV, 
manually move the 5 normal CSVs into the Normal_CSV folder, and then move
the 5 tumor CSVs into the Tumor_CSV folder, this can be done 
by using the search bar in the Finder(Mac) or Folder(Windows) app. 
The script can then point to the directory (similar to HW1) to read and merge the files within, 
using a function within the pandas (pd) package.

Reading in a CSV file Example:
DataFrame1 = pd.read_csv("DataFrame1.csv")

Merging Example:
newDataFrame = pd.concat(DataFrame1, DataFrame2, axis=0) 
'''
normal1 = pd.read_csv('Homework_3/Normal_CSV/ab76efd7-0859-4bb5-8da7-a2185ffc0567_normal.csv')
normal2 = pd.read_csv('Homework_3/Normal_CSV/ab6504e6-37e4-451a-9530-f9aa88a18263_normal.csv')
normal3 = pd.read_csv("Homework_3/Normal_CSV/afcba237-47af-42fa-b624-f664773abdef_normal.csv")
normal4 = pd.read_csv("Homework_3/Normal_CSV/b1dbbd1e-f48a-4bcc-9618-6c89c5c98f51_normal.csv")
normal5 = pd.read_csv("Homework_3/Normal_CSV/b2a8da4b-6c32-4afb-a23d-bd14f858be58_normal.csv")

normal_df = pd.concat([normal1, normal2, normal3, normal4, normal5], axis = 0)


tumor1 = pd.read_csv('Homework_3/Tumor_CSV/ab76efd7-0859-4bb5-8da7-a2185ffc0567_tumor.csv')
tumor2 = pd.read_csv('Homework_3/Tumor_CSV/ab6504e6-37e4-451a-9530-f9aa88a18263_tumor.csv')
tumor3 = pd.read_csv('Homework_3/Tumor_CSV/afcba237-47af-42fa-b624-f664773abdef_tumor.csv')
tumor4 = pd.read_csv('Homework_3/Tumor_CSV/b1dbbd1e-f48a-4bcc-9618-6c89c5c98f51_tumor.csv')
tumor5 = pd.read_csv('Homework_3/Tumor_CSV/b2a8da4b-6c32-4afb-a23d-bd14f858be58_tumor.csv')

tumor_df = pd.concat([tumor1, tumor2, tumor3, tumor4, tumor5], axis = 0)



# Function adds in alt_seq column to, input is a dataframe and function returns a dataframe
def addALT_Seq(csv):
    alt = []
    for row in range(csv.shape[0]):
        ref_seq = csv["ref_seq"].iloc[row]
        if ref_seq == csv["var_seq1"].iloc[row]:
            alt.append(csv["var_seq2"].iloc[row])
        else:
            alt.append(csv["var_seq1"].iloc[row])
    csv.insert(csv.shape[1], "alt_seq", alt)
    return csv

'''
INSERT CODE HERE
Homework 3.A
(ii) Using the output from A(i), run the addALT_Seq() function:
Example:
newDataFrame_withALTseq = addALT_Seq(NewDataFrame)
'''
normal_withALTseq = addALT_Seq(normal_df)
tumor_withALTseq = addALT_Seq(tumor_df)

'''
INSERT CODE HERE
Homework 3.A
(iii) Using the output from A(ii), remove duplicates based on the given columns:
[“chrom”, “left”, “ref_seq”, “alt_seq”, “Patient_ID”]
Save the two DataFrames as: Final_Normal and Final_Tumor

Remove Duplicates Example:
Final = newDataFrame_withALTseq.drop_duplicates(columns)
'''
columns = ["chrom", "left", "ref_seq", "alt_seq", "Patient_ID"]
Final_Normal = normal_withALTseq.drop_duplicates(subset=columns)
Final_Tumor = tumor_withALTseq.drop_duplicates(subset=columns)

normal_csv = Final_Normal.to_csv('/Users/yeseniajuarez/Documents/genomic_analysis/Homework_3/normal_merged.csv')
tumor_csv = Final_Tumor.to_csv('/Users/yeseniajuarez/Documents/genomic_analysis/Homework_3/tumor_merged.csv')
'''
OUTPUT CHECK
Homework 3.A
(iv) Run the lines below:
'''
print("The number of (Rows, Columns) in Tumor:")
print(Final_Tumor.shape)
print("The number of (Rows, Columns) in Normal:")
print(Final_Normal.shape)

# 1. filter SNPNexus tsv files
tumor_file_pervariant = pd.read_csv("Homework_4/HW4_tumor_pervariant.tsv", sep='\t')
normal_file_pervariant = pd.read_csv("Homework_4/HW4_normal_pervariant.tsv", sep='\t')

tumor_filtered = tumor_file_pervariant[~tumor_file_pervariant['id'].str.endswith(':-1')]
normal_filtered = normal_file_pervariant[~normal_file_pervariant['id'].str.endswith(':-1')]

tumor_filtered.to_csv("Homework_4/HW4_tumor_pervariant_filtered.tsv", sep='\t', index=False)
normal_filtered.to_csv("Homework_4/HW4_normal_pervariant_filtered.tsv", sep='\t', index=False)

print(tumor_file_pervariant.shape)
print(tumor_filtered.shape)

# 2. VCF file prep for FATHMM, drop downregulated -1

tumor_vcf = pd.read_csv("/Users/yeseniajuarez/Documents/genomic_analysis/Homework_4/gen_coords_85db7982.vcf", sep='\t', skiprows=8)
normal_vcf = pd.read_csv("/Users/yeseniajuarez/Documents/genomic_analysis/Homework_4/gen_coords_f9fb2794.vcf", sep='\t', skiprows=8)

tumor_vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
normal_vcf.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

# filter to remove downrregulated -1
tumor_vcf_filtered = tumor_vcf[~tumor_vcf['ID'].str.endswith(':-1')]
normal_vcf_filtered = normal_vcf[~normal_vcf['ID'].str.endswith(':-1')]

tumor_vcf_filtered.to_csv("Homework_4/tumor_vcf_filtered.vcf", sep='\t', index=False)
normal_vcf_filtered.to_csv("Homework_4/normal_vcf_filtered.vcf", sep='\t', index=False)

# read FATHMM files
fathmm_tumor = pd.read_csv("Homework_4/tumor_fathmm.txt", sep='\t')
fathmm_normal = pd.read_csv("Homework_4/normal_vcf.txt", sep='\t')

# reset index
Final_Tumor_reset = Final_Tumor.reset_index(drop=True)
tumor_filtered_reset = tumor_filtered.reset_index(drop=True)
fathmm_tumor_reset = fathmm_tumor.reset_index(drop=True)

Final_Normal_reset = Final_Normal.reset_index(drop=True)
normal_filtered_reset = normal_filtered.reset_index(drop=True)
fathmm_normal_reset = fathmm_normal.reset_index(drop=True)

# concatenate all 3 data frames
tumor_final = pd.concat([Final_Tumor_reset, tumor_filtered_reset, fathmm_tumor_reset], axis=1)
normal_final = pd.concat([Final_Normal_reset, normal_filtered_reset, fathmm_normal_reset], axis=1)

# Save final files
tumor_final.to_csv("Homework_4/tumor_final.csv", index=False)
normal_final.to_csv("Homework_4/normal_final.csv", index=False)

print(f"Final tumor file: {tumor_final.shape}")
print(f"Final normal file: {normal_final.shape}")



