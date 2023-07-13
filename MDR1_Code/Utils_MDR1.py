import pandas as pd
import numpy as np
import RNA
import re
import os

variant_info = {
    1: {
        "variant_name": "T1236C",
        "aa_position": 412, 
        "cds_position": 1236,
        "genome_position": 87550285,
        "codon_before": "GGT",
        "codon_after": "GGC",
        "change_from": "T",
        "change_to": "C"},
    2: {
        "variant_name": "T2677G",
        "aa_position": 893,
        "cds_position": 2677,
        "genome_position": 87541302,
        "codon_before": "TCT",
        "codon_after": "GCT",
        "change_from": "T",
        "change_to": "G"},
    3: {
        "variant_name": "T3435C",
        "aa_position": 1145,
        "cds_position": 3435,
        "genome_position": 87509329,
        "codon_before": "ATT",
        "codon_after": "ATC",
        "change_from": "T",
        "change_to": "C"}}


SynonymousCodons = {
  "C": ["TGT", "TGC"],
  "D": ["GAT", "GAC"],
  "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
  "Q": ["CAA", "CAG"],
  "M": ["ATG"],
  "N": ["AAC", "AAT"],
  "P": ["CCT", "CCG", "CCA", "CCC"],
  "K": ["AAG", "AAA"],
  "*": ["TAG", "TGA", "TAA"],
  "T": ["ACC", "ACA", "ACG", "ACT"],
  "F": ["TTT", "TTC"],
  "A": ["GCA", "GCC", "GCG", "GCT"],
  "G": ["GGT", "GGG", "GGA", "GGC"],
  "I": ["ATC", "ATA", "ATT"],
  "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
  "H": ["CAT", "CAC"],
  "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
  "W": ["TGG"],
  "V": ["GTA", "GTC", "GTG", "GTT"],
  "E": ["GAG", "GAA"],
  "Y": ["TAT", "TAC"],"-": ["---"]}

reverse_dict = {value: key for key in SynonymousCodons for value in SynonymousCodons[key]}
 
nucs_dict = {"A": 0, "C": 1, "G": 2, "T": 3} #concensus of positions for the rest of the code: A will always be the 0 positon in the matrices and so on

''' 
We use tAI values of *healthy tissues* but their names in our dictionary of tAI weights is according to the TCGA project in the current tissue. 
This dictionary maps between the name of the TCGA project to the tissue 
'''
tissues_titles_dict = {
    "KIRC": "Kidney",
    "LIHC": "Liver",
    "GBM_": "Brain",
    "COAD": "Colon"
}
     

'''This function returns a (4,3) matrix holding all possible synonymous substitutions of a codon. Each row of the
matrix corresponds to a single nucleotide in the ACGT alphabet (in this alphabetical order) and each column
correseponds to a nucloetide position in the codon (first, middle, last).  
For example, TGT and TGC are the only codons which code for Cys(C). So, for "TGT" we will get this matrix:
[0,0,0],
[0,0,1],
[0,0,0],
[0,0,0]

This means that the only possible synonymous substitution for this codon is the third codon changing to a "C". ''' 

def possible_syn_muts_per_codon(cur_codon:str, AA_to_codon_dict:dict, nucs_to_position_dict:dict) -> np.ndarray:
    mat = np.zeros((4, 3)) 
    for i in range(3): #iterate over positions in codons
        for j in nucs_to_position_dict.keys():
            if cur_codon[i] != j: #if they are equal than we are not making any change
                changed_codon = cur_codon
                temp = list(changed_codon) #strings are emutable
                temp[i] = j #mutate
                changed_codon = ''.join(temp)
                if AA_to_codon_dict[cur_codon] == AA_to_codon_dict[changed_codon]: #check if the mutation was synonymous
                    mat[nucs_to_position_dict[j],i] = 1
    return(mat)

'''Creating a dictionary that holds the synonymous substitution matrices for all codons (general, not
specific to any gene)'''

codons_syn_maps_dict = {}
for codon in reverse_dict.keys():
    if codon != "---":
        mat_of_codon = possible_syn_muts_per_codon(codon, reverse_dict, nucs_dict)
        codons_syn_maps_dict[codon] = mat_of_codon
    
    
'''This function creates a matrix of shape (4,len_cds) which holds the possible synonymous substitutions of our
specific cds'''
def get_possib_syn_sub_per_positon(nt_CDS:str ,codons_syn_maps_dict:dict) -> np.ndarray:

    # split to codons and replace every codon with its matrix of possible synonymous substitutions
    split_to_codons = [nt_CDS[i:i+3] for i in range(0, len(nt_CDS), 3)] #split cds string to codons
    split_to_codons = np.asarray(split_to_codons,dtype = object)
    for codon, syn_mat in codons_syn_maps_dict.items():
        indices = np.where(split_to_codons == codon)[0] #all locations of this codon in the CDS of the current gene 
        for index in indices:
            split_to_codons[index] = syn_mat #change from the codon to its synonymous substitution matrix
    possible_syn_replacements_for_gene = np.concatenate(split_to_codons,axis = 1) #concatinate all matrices together
    #now we got a matrix of shape "len_cds" that contains the possible synonymous substitutions of each position. 
    return(possible_syn_replacements_for_gene)
    


'''Perform a single mutation in a sequence'''
def mutate_cds_sequence(sequence, position, change_to):
    position = position - 1 #change to 0-based
    temp = list(sequence)
    temp[position] = change_to
    mutant_seq = "".join(temp)
    return(mutant_seq)

'''This function calculates the mfe of a single sequence using the Vienna package'''
def calc_windows_mfe(nt_CDS, window_size):
    num_windows = len(nt_CDS) - window_size + 1
    mfe_windows = np.empty((1,num_windows))
    mfe_windows[:] = np.nan

    for window in range(num_windows):
        local_seq = nt_CDS[window: window + window_size]
        assert(len(local_seq) == 39)
        (_, mfe) = RNA.fold(local_seq)
        mfe_windows[0,window] = mfe
    return(mfe_windows)

'''We need to sum the columns of "mfe_windows_padded" and its shifted versions. We need for
nan + nan to be equal nan but the numpy function np.nansum returns zero in that case. 
(This is bad for us because a 0 score means positions with weak folding and this is not true for
nan positions). So, we create a wrapper for that function, that returns a nan if all the values
summed are equal nan.'''

def nansumwrapper(vectors_to_sum, axis):
    nan_positions = np.isnan(vectors_to_sum).all(axis) #positions where both vectors that are being summed are equal nan. 
    res = np.nansum(vectors_to_sum, axis) #the result of the summation, there are zeros where only nans are summed
    res[nan_positions] = np.nan #insert nans in those locations
    return res

'''In this function we perform the transition from mfe of windows to mfe of positions in the cds. The mfe score of a position 
is the mean of the scores of all windows the position is in. We are creating "window_size - 1" shifted versions of the mfe_window
results vector and so we take the mean of the windows in a vectorized manner. 

example: 
let's say we have an original sequence of length 10, and the window size is 4. Then the number of windows is 7
(length_cds - window_size + 1). 
mfe_windows = [-1,-3,-2,-2,-4,-3,-2]
according to the rule, the results should be:
mfe_position_0 = mean(mfe_windows[0]) = -1 (position 0 is only in window 0)
mfe_position_1 = mean(mfe_windows[0,1]) = ((-1-3)/2) = -2
mfe_position_2 = mean(mfe_windows[0,1,2]) = ((-1-3-2)/3) = -2
mfe_position_3 = mean(mfe_windows[0,1,2,3]) = ((-1-3-2-2)/4) = -2
mfe_position_4 = mean(mfe_windows[1,2,3,4]) = ((-3-2-2-4)/4) = -2.75
mfe_position_5 = mean((mfe_windows[2,3,4,5])) = ((-2-2-4-3)/4) = -2.75
mfe_position_6 = mean((mfe_windows[3,4,5,6])) = ((-2-4-3-2)/4) = -2.75
mfe_position_7 = mean((mfe_windows[4,5,6])) = ((-4,-3,-2)/3) = -3
mfe_position_8 = mean((mfe_windows[5,6])) = ((-3,-2)/2) = -2.5
mfe_position_9 = mean((mfe_windows[6])) = -2

let's show that we are getting the same result using our method (but fast)
pad the mfe_windows vector: add nans to the end to make it in the length of "num_positions"
mfe_windows = [-1,-3,-2,-2,-4,-3,-2] (length 7) -> mfe_windows_padded = [-1,-3,-2,-2,-4,-3,-2,nan,nan,nan] (length 10)
shifted_1 = [nan,-1,-3,-2,-2,-4,-3,-2,nan,nan] (one score was moved from the end to the beginning)
shifted_2 = [nan,nan,-1,-3,-2,-2,-4,-3,-2,nan] (two scores were moved from the end to the beginning)
shifted_3 = [nan,nan,nan,-1,-3,-2,-2,-4,-3,-2] 
how many shifted versions do we need? window_size -1

now, we summing mfe_windows,shifted_1,shifted_2 and shifted_3. 
(in the code we do this one after the other so we dont have to keep all the matrices, because this is memory inefficient)
We take the average of each column:
mean_col0 = mean([-1,nan,nan,nan]) = -1
mean_col1 = mean([-3,-1,nan,nan]) = -2
mean_col2 = mean([-2,-3,-1,nan]) = -2
mean_col3 = mean([-2,-2,-3,-1]) = -2
mean_col4 = mean([-4,-2,-2,-3]) = -2.75
mean_col5 = mean([-3,-4,-2,-2]) = -2.75
mean_col6 = mean([-2,-3,-4,-2]) = -2.75
mean_col7 = mean([nan,-2,-3,-4]) = -3
mean_col8 = mean([nan,nan,-2,-3]) = -2.5
mean_col9 = mean([nan,nan,nan,-2]) = -2
(Notice that we are summing over the not-nan positions so we have to keep track of the amount of nans per position). 

we can see that we got the same results! but the second way can be vectorized so we will use it.
'''
def calc_mfe_per_position(mfe_windows:np.ndarray, data:str, padding_amount:int, window_size:int) -> np.ndarray:
    #perform the padding
    if data == 'original':
        mfe_windows_padded = np.pad(mfe_windows, ((0,0),(0, padding_amount)), 'constant', constant_values=(np.nan))
    else:
        mfe_windows_padded = np.pad(mfe_windows, ((0,0),(0,0),(0, padding_amount)), 'constant', constant_values=(np.nan))
    #perform the averaging. #each time, we shift and sum. later we devide by the number of windows the position was in. 
    mfe_per_position = mfe_windows_padded.copy() #initialize
    nans_per_position = np.isnan(mfe_windows_padded).astype(int) #keep track of nans in order to devide later by num of non-nan positions
    axis = 1 if data == 'original' else 2
    for i in range(1,window_size):
        shifted = np.roll(mfe_windows_padded, i, axis = axis) #create the next shifted version
        nans_per_position = nans_per_position + np.isnan(shifted).astype(int) #where are the nans in this shifted version? 
        mfe_per_position = nansumwrapper([mfe_per_position, shifted], axis = 0)#sum the prior columns of the prior result with this shifted version

    not_nans_per_position = window_size - nans_per_position #for each column, how many elements were not equal to nan?
    mfe_per_position = mfe_per_position/not_nans_per_position #for each column, devide by the number of not nan elements
    return(mfe_per_position)

''' Related to the fptc '''
def get_codon_and_freq(freq_data_per_codon):
    
    codon_and_freuqency = re.search("[A-Z](.*?)\(", freq_data_per_codon)
    codon_and_freuqency = codon_and_freuqency[0][:-1].split(" ")
    codons = codon_and_freuqency[0]
    freq1000 = codon_and_freuqency[-1]
    
    return(codons, freq1000)

''' This function calculates and empirical p-value given the delta mfe caused by the true mutation and 
a vector of the delta mfes caused by the random mutations '''
def get_pvalue(true_delta: float, random_deltas: np.array) -> float:
    
    pval = 1 - ((np.sum(true_delta > random_deltas)) / len(random_deltas)) 
    return(pval)

''' The maf files have varying number of header rows. This function gets the number of header lines of a specific file, 
so we can input it to pd.read_csv and get the maf file as a table'''
def get_rows_to_skip(path):
    header_counter = 0
    with open(path, "r") as a_file:
        for line in a_file:
            if line[0]!= '#': #header rows start with the "#" sign
                break
            header_counter += 1
    return(header_counter)


''' Get all TCGA variant data for a specific gene'''
def get_cancer_muts_cur_gene(gene_name):
    cancerous_muts_path = f"../../../../../tamir1/cancer_proj/gdc_db/data/filtered_feb_2021/AllGenes/{gene_name}/"
    #iterate over all the MAF files in the folder of this gene and get all synonymous mutations
    cancer_muts_cur_gene = pd.DataFrame()
    for filename in os.listdir(cancerous_muts_path):
        if filename.endswith(".maf"):
            path_cur_cancer_type = f"{cancerous_muts_path}{filename}"
            num_rows_to_skip = get_rows_to_skip(path_cur_cancer_type) #the number of header notes is different for each file
            temp_df = pd.read_csv(path_cur_cancer_type , sep = "\t", skiprows = num_rows_to_skip)
            #temp_df = temp_df[temp_df['Variant_Classification'] == 'Silent'] #keep only synonymous variants
            temp_df['Cancer_Type'] = filename.split(".")[0] #keep the data of the cancer type
            cancer_muts_cur_gene = pd.concat([cancer_muts_cur_gene, temp_df], sort=False)
    return(cancer_muts_cur_gene)

def get_muts_single_patient(case_id, patient_cancer_type):

    cur_patient_muts = pd.DataFrame()

    root_path = "../../../../../tamir1/cancer_proj/gdc_db/data/filtered_feb_2021/AllGenes/"
    for root, dirs, files in os.walk(root_path):
        for name in files:
            if patient_cancer_type in name:
                path_to_df = os.path.join(root, name)
                num_rows_to_skip = get_rows_to_skip(path_to_df) #the number of header notes is different for each file
                temp_df = pd.read_csv(path_to_df , sep = "\t", skiprows = num_rows_to_skip)
                temp_df = temp_df[temp_df["case_id"] == case_id]
                cur_patient_muts = pd.concat([cur_patient_muts, temp_df], sort=False)
                
# def get_num_muts_single_patient(case_id, patient_cancer_type):

#     counter = 0
    
#     root_path = "../../../../../tamir1/cancer_proj/gdc_db/data/filtered_feb_2021/AllGenes/"
#     for root, dirs, files in os.walk(root_path):
#         for name in files:
#             if patient_cancer_type in name:
#                 path_to_df = os.path.join(root, name)
#                 num_rows_to_skip = get_rows_to_skip(path_to_df) #the number of header notes is different for each file
#                 temp_df = pd.read_csv(path_to_df , sep = "\t", skiprows = num_rows_to_skip)
#                 temp_df = temp_df[temp_df["case_id"] == case_id]
#                 counter += temp_df.shape[0]

    return(counter)

def reverse_complement(s: str, complement: dict = None) -> str:
    """Performs reverse-complement of a sequence. Default is a DNA sequence."""
    if complement is None:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    s_rev = s[::-1]
    lower = [b.islower() for b in list(s_rev)]
    bases = [complement.get(base, base) for base in list(s_rev.upper())]
    rev_compl = ''.join([b.lower() if l else b for l, b in zip(lower, bases)])
    return rev_compl


