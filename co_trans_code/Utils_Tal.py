#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import numpy as np
import RNA
import pandas as pd
import re
'''
Here we define two dictionaries and a list of valid AAs:
SynonymousCodons - a dictionary whose keys are amino-acids and whose values are codons.
reverse_dict - a dictionary whose keys are codons and whose values are amino-acids. 
'''


AAs = ["C","D","S","Q","M","N","P","K","T","F","A","G","I","L","H","R","W","V","E","Y"]

#the non-silent categories in TCGA
type_NonSilent = ['Missense_Mutation', 'In_Frame_Del', 'Nonsense_Mutation', 'In_Frame_Ins', 'Nonstop_Mutation','Frame_Shift_Del','Frame_Shift_Ins', 'Translation_Start_Site']

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
nucs_dict_reverse = {0:"A", 1:"C", 2:"G", 3:"T"}

''' This function recieves an array of aa positions and turns it to an array of nt positions.'''
def aa_positions_to_nt_positions(aa_positions):
    if np.all(np.isnan(aa_positions)):
        return(aa_positions)
    else:
        nt_positions = []
        first_nt = aa_positions * 3
        second_nt = first_nt + 1
        third_nt = second_nt + 1
        nt_positions.extend(first_nt)
        nt_positions.extend(second_nt)
        nt_positions.extend(third_nt)
        nt_positions = np.array(nt_positions)
        return(nt_positions)

''' This function recieves a single nucleotide position and turns it to an amino-acid position'''
def nt_position_to_aa(nt_position):
    if np.all(np.isnan(nt_position)):
        return(nt_position)
    else:
        aa_position = int(np.floor(nt_position / 3))
        return(aa_position)
    
    
''' Given a nucleotide position, get the respective codon position and the position of the nt in the codon itself ''' 
def get_position_in_codon_and_codon_position(nt_position):
    aa_position = nt_position_to_aa(nt_position)
    pos_relative_codon = nt_position % (aa_position * 3)
    return(aa_position, pos_relative_codon)
    

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
    
'''A dictionary that holds the synonymous substitution matrices for all codons.'''
codons_syn_maps_dict = {}
for codon in reverse_dict.keys():
    if codon != '---':
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
    

'''This function takes a single randomization as if it was the original gene and compares it to the rest
of the randomizations. It checks how many positions are chosen as significant (extremly low/high cai)
according to the "significance_threshold" chosen. This is the estimated number of false positive results 
that we expect to get when comparing to the original gene. The function return the number of false
positive results expected for high cai and low cai'''
def evaluate_false_positive_rate(score_rands:np.ndarray, significance_threshold:float) -> [int, int]:
    
    score_single_rand, score_rest_rands = score_rands[0], score_rands[1:] #split a randomizations from the batch
    num_rest_rands = score_rest_rands.shape[0]

    single_is_larger = np.sum(score_single_rand > score_rest_rands , axis = 0) / num_rest_rands #the single randomization has a higher CAI than how many randomizations?
    num_fp_larger = np.sum(single_is_larger >= significance_threshold)
    single_is_smaller = np.sum(score_single_rand < score_rest_rands , axis = 0) / num_rest_rands
    num_fp_smaller = np.sum(single_is_smaller >= significance_threshold)
    
    return(num_fp_larger, num_fp_smaller)    

'''Gets 3 lists (groups of genes), calculates all the needed groups to plot a venn diagram of 3 groups (gene in a but not in b, c, etc.)
plots the diagram'''
def venn_diagram(a, b, c, labels, title):
    a_not_bc = len([gene for gene in a if gene not in b and gene not in c])
    b_not_ac = len([gene for gene in b if gene not in a and gene not in c])
    c_not_ab = len([gene for gene in c if gene not in a and gene not in b])
    ab_not_c = len([gene for gene in a if gene in b and gene not in c])
    ac_not_b = len([gene for gene in a if gene in c and gene not in b])
    bc_not_a = len([gene for gene in b if gene in c and gene not in a])
    abc = len([gene for gene in a if gene in b and gene in c])
    plt.title(title)
    venn3(subsets = (a_not_bc, b_not_ac, ab_not_c, c_not_ab, ac_not_b, bc_not_a, abc), set_labels = labels, alpha = 0.5);   
    
    
'''Perform a single mutation in a sequence
* sequence: the original sequence
* position: position to mutate
* change_to: change to this nucleotide (either A, C, G or T)
* coordiante_system: the position is given in a 0-based or 1-based meathod'''
def mutate_cds_sequence(sequence, position, change_to, coordiante_system):
    if coordiante_system == 1:
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



'''Get the human codon frequencies per 1000 and input to dictionary'''


''' The downloaded data is in a form that is not tabular and not convinent for programming. Here we simply change the from. '''

def get_codon_and_freq(freq_data_per_codon):
    
    codon_and_freuqency = re.search("[A-Z](.*?)\(", freq_data_per_codon)
    codon_and_freuqency = codon_and_freuqency[0][:-1].split(" ")
    codons = codon_and_freuqency[0]
    freq1000 = codon_and_freuqency[-1]
    
    return(codons, freq1000)



codon_usage_dict = pd.read_csv("../Data/Human_codon_frequency.txt", skiprows = 5, header = None) #human
codons_and_freq1000 =  codon_usage_dict[0].apply(lambda x: get_codon_and_freq(x))
dict_freq1000 = dict((codon, freq) for codon, freq in codons_and_freq1000)
