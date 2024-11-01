# pylint: disable=line-too-long,invalid-name,pointless-string-statement,too-many-arguments

"""
Enformer Utils

Credits to Yoram Zarai from the Tuller lab at Tel Aviv university who wrote almost all of the code in this file.   

"""
import sys
import pathlib
from typing import List, Dict, Optional, Callable
import pickle
import matplotlib.pyplot as plt
from matplotlib import ticker
import tensorflow_hub as hub
import numpy as np
import seaborn as sns
import pandas as pd

# make sure PYTHONPATH is set, containing the Utils/ folder.
from Fasta_segment import Fasta_segment

#import translation_py36 as tran

# due to the old python, can not use instead Utils/mysys_init.py
if  sys.platform == 'linux':
    MACHINE = 'POWER'
    ROOT_path = '../Projects/MDR1/'
    chromosome_path = pathlib.Path('../Data/Genomes/Human/human_hg38/Chromosome')
else:
    raise Exception(f"{__name__}: Platform {sys.platform} is not supported !!")


# Enformer constants
# ==================
Enformer_model_url: str = "https://tfhub.dev/deepmind/enformer/1"
Enformer_input_len: int = 393_216  # actual input size (including flanking NTs)
Enformer_predict_len: int = 114_688  # the output prediction corresponds to the center 114_688 of the 393_216 input
Enformer_num_tracks: int = 5313
Enformer_out_len: int = 896  # = 114688/128 referred to here as bins. This is the number of outputs (bins) per track
Enformer_NTs_per_bin = int(Enformer_predict_len / Enformer_out_len)  # = 128. Each bin output corresponds to 128 NTs

# Enformer related data
Enformer_target_info_file = '../Data/enformer_targets_df.pkl'  # dataframe containing tracks information (generated by enformer_tensorflow_usage.ipynb)


# Sequence Functions
# ==================
def pull_fasta_sequence(fasta_file: str, start_loc: int, end_loc: int) -> str:
    """
    Given a start and end 1-based positions, the functions generates the corresponding sequence from the
    chromosome fasta file.

    start_loc - 1-based index of first NT of the sequence of interest.
                1 corresponds to the first NT in the fasta file.
    end_loc - 1-based index of last NT of the sequence of interest.
    """
    return Fasta_segment().read_segment(fasta_file, start_loc-1, end_loc-start_loc+1)

def get_chrm_file_name_mac(chrm_num: str) -> str:
    """
    Chromosome fasta file name on my MAC.
    """
    return "chr" + chrm_num + ".fa"


def get_chrm_file_name_power(chrm_num: str) -> str:
    """
    Chromosome fasta file name on Power.
    """
    return "Homo_sapiens.GRCh38.dna_sm.chromosome." + chrm_num + ".fa"


get_chrm_file_name: Callable[[str], str] = get_chrm_file_name_mac if MACHINE == 'MAC' else get_chrm_file_name_power


def generate_mut_variant(seq: str, var_type: str,
                         start_pos: int, end_pos: int,
                         mut: str) -> str:
    """
    This function returns a mutated sequence based on the input
    sequence and mutation information.

    start_pos: start position of where to mutate the sequence seq. 0-base number,
               i.e. start_pos=0 implies the first character is seq.
    end_pos: end position of where to mutate the sequence seq. 0-base number.

    var_type: 'SNP', 'DNP', 'TNP', 'ONP', 'DEL' or 'INS'.

    mut: the SNP/DNP/TNP/ONP mutant (in case of var_type in ['SNP', 'DNP', 'TNP', 'ONP']),
    or the sequence to insert in case of 'INS' var_type. mut is not used in case of var_type=='DEL'.
    """
    if var_type in ['SNP', 'DNP', 'TNP', 'ONP']:
        '''
        In case of SNP/DNP/TNP/ONP, all characters from start_pos to end_pos are altered.
        Thus, both start_pos and end_pos in seq are altered.
        '''
        return seq[:start_pos]+mut+seq[end_pos+1:]
    elif var_type == 'INS':
        '''
        In case of INS, the inserted mut is between start_pos and end_pos.
        Thus, both start_pos and end_pos positions in seq are not altered.
        '''
        return seq[:start_pos+1]+mut+seq[end_pos:]
    elif var_type == 'DEL':
        '''
        In case of DEL, all characters from start_pos to end_pos are removed.
        '''
        return seq[:start_pos]+seq[end_pos+1:]
    else:
        print(f"generate_mut_variant: var_type {var_type} not supported !!")
        return ''


# Enformer Functions
# ==================
def load_enformer():
    """loads the trained model"""
    return hub.load(Enformer_model_url).model


def enformer_predict(enformer, one_hot: np.ndarray) -> dict:
    """Given a one-hot encoded sequence (of shape (Enformer_input_len, 4)), runs Enformer."""
    assert one_hot.shape == (Enformer_input_len, 4), f"Enformer one-hot encoded input sequence shape must be ({Enformer_input_len}, 4), but is {one_hot.shape} !!"
    #print("DEBUD CODE in enfoemr_predict() in enformer_utils.py !!!!!!!!!!!!")
    #return {'a': 1, 'b': 2, 'human': [1,2,3,4,5]}
    return enformer.predict_on_batch(one_hot[np.newaxis])


transform_CAGE_output: Callable[[np.ndarray], np.ndarray] = lambda x: np.log10(1.0 + x)


def interval_to_vals(interval: str) -> dict:
    """Given 'chrm:start-end' string, returns chrm, int(start), and int(end)."""
    chrm, region = interval.split(':')
    start, end = map(int, region.split('-'))
    return {'chrm': chrm, 'start': start, 'end': end}


def get_ext_interval(interval: str, ex_size: int) -> str:
    """
    Converts the predicted interval (i.e. 'interval') string to the extended input interval
    (i.e. the actual input range to the Enformer) string.
    For interval=chrm:x-y, then y-x+1 = 116888, so x [y] is the
    first [last] NT of the prediction range.

    If c is the center position (i.e. c = (start+end)//2), then
    1. The predicted interval (i.e. the 'interval' input) corresponds to:
        <57,343 NTs><c><57,344 NTs>, i.e. 114,688 NTs.
    2. The corresponding extended input interval (i.e. the actual Enformer input) is:
        <196,607 NTs><c><196,608 NTs>, i.e. 393,216 NTs.

    The code below is general for any even or odd predicted size and extended size.
    """
    intr = interval_to_vals(interval)
    c = (intr['start'] + intr['end']) // 2
    # the last NT in the first half is the "center", as the interval (and the
    # extended interval) size is even. Thus adding "+1" to the first half.
    nstart, nend = c - (ex_size // 2) + 1 - ex_size % 2, c + (ex_size // 2)
    return f"{intr['chrm']}:{nstart}-{nend}"


def enformer_intervals(mid_position: int, chrm_str: str) -> tuple:
    """
    Computes the predicted interval and the input interval.
    chrm_str: the chromosome number, e.g. 'chr7'.

    The mid_position is the position in the "middle" of the sequence to be predicted (predicted interval), 
    i.e. the predicted sequence contains 57,343 upstream NTs, mid_position NT, and 57,344 downstream NTs. 

    returns predict_interval_start, predict_interval_end, predict_interval_str,
    input_interval_start, input_interval_end, input_interval_str
    """
    # predicted range
    predt_start = mid_position - (Enformer_predict_len // 2) + 1
    predt_end = predt_start + Enformer_predict_len - 1 
    predt_interval = f"{chrm_str}:{predt_start}-{predt_end}"
    
    # input range
    input_interval = get_ext_interval(predt_interval, Enformer_input_len)
    vals = interval_to_vals(input_interval)

    return predt_start, predt_end, predt_interval, vals['start'], vals['end'], input_interval


def plot_tracks(tracks: dict, axes, interval: str,
                rev: bool = False,
                xtick_format: str = ',.0f',
                xticks_num: int = 8):
    """
    Prints track per sub-plot.
    tracks: a dict with key contains a string dscription, and value of numpy vector.

    if rev==True, we flip the plot, i.e. we show the prediction relative to the
    positive strand.
    """
    d_interval = interval_to_vals(interval)
    # if rev==True, the input to enformer is the reveresed-complement of the positive strand.
    # The outputs (tracks) correspond to the reversed-complement. We flip it to show
    # the prediction relative to the positive strand.
    srt, end = (d_interval['end'], d_interval['start']) if rev else (d_interval['start'], d_interval['end'])
    predct_rng = np.abs(end - srt) + 1 #d_interval['end'] - d_interval['start']
    for ax, (title, y) in zip(axes, tracks.items()):
        # this flips the plot in case rev==True
        ax.fill_between(np.linspace(srt, end, num=len(y)), y)
        ax.set_title(title)
        sns.despine(top=True, right=True, bottom=True)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(predct_rng/(xticks_num-1)))
        ax.xaxis.set_major_formatter(lambda x, pos: f'{x:{xtick_format}}')
        ax.set_xlabel(interval)
    plt.tight_layout()


def pos2bin(start_pos: int, pos: int, rev: bool = False) -> tuple:
    """
    Given a coordinate position 'pos' in the predicted
    positive strand input range of Enformer, the function returns
    the corresponding Enformer's output 1-based bin number, and the 1-based NT offset within the bin.

    A bin is the 1-based index of the Enformer output, where each bin corresponds to 128 consecutive NTs.
    There is a total of 896 output bins.

    Note that the bins correspond to the reveresed-complement sequence in case of rev == True.

    start_offset is the coordinate of the first NT in the 114,688 positive strand predicted range of the Enformer.
    This implies that start_pos <= pos.
    """
    assert pos >= start_pos, f"{pos=:,} must be greater or equal to {start_pos=:,} !!"
    assert (pos-start_pos) < Enformer_predict_len, f"{pos=:,} is outside the predicted range of Enformer (given {start_pos=:,}) !!"

    offset = pos - start_pos
    # 1-based, relative to the positive strand
    bin_num, bin_offset = (offset // Enformer_NTs_per_bin) + 1, (offset % Enformer_NTs_per_bin) + 1
    return (Enformer_out_len-bin_num+1, Enformer_NTs_per_bin-bin_offset+1) if rev else (bin_num, bin_offset)


def bin2range(start_pos: int, bin_num: int, rev: Optional[bool] = False) -> tuple:
    """
    Given the 1-based bin number bin_num, the function returns
    the corresponding coordinate range that the bin covers.

    A bin is the 1-based index of the Enformer output, where each bin corresponds to 128 consecutive NTs.
    There is a total of 896 output bins.

    Note that if rev==True, 'bin_num' is the bin number in the reveresed-complement strand (i.e. the negative strand)

    start_offset is the coordinate of the first NT in the 114,688 positive-strand predicted range of the Enformer.
    """
    # the bin number in the positive strand sequence
    assert 1 <= bin_num <= Enformer_out_len, f"{bin_num=} must be >=1 and <= {Enformer_out_len} !!"

    pos_strand_bin_num = (Enformer_out_len-bin_num+1) if rev else bin_num

    # the range corresponds to coordinates in the positive strand, i.e. start < end
    start, end = start_pos + (pos_strand_bin_num-1)*Enformer_NTs_per_bin, start_pos + pos_strand_bin_num*Enformer_NTs_per_bin - 1

    return (end, start) if rev else (start, end)


def adjust_seq_size(seq: str, end_pos: int, chrm_path: str, req_seq_len: int, rev: bool) -> str:
    """
    Trims input sequence if too large, or adds chromosome sequence if too short. Used in case
    of INS and DEL mutations.

    rev == True implies that the sequence was reveresed complemented. Recall that the
    mutation is added to the positive strand sequence.
    """
    o_seq = seq
    if len(o_seq) > req_seq_len:
        # INS mutation
        # INS mutation added to the positive strand. Thus, if rev=True, need to trim the beginning
        o_seq = o_seq[len(o_seq)-req_seq_len:] if rev else o_seq[:req_seq_len]
    elif len(o_seq) < req_seq_len:
        # DEL mutation
        extra = req_seq_len - len(o_seq)
        extra_seq = pull_fasta_sequence(chrm_path, end_pos + 1, end_pos + extra).upper()
        o_seq = reverse_complement(extra_seq) + o_seq if rev else o_seq + extra_seq
    return o_seq


def enformer_predict_both_strands(ref_seq: str, var_type: str, rel_sp: int, rel_ep: int, mut_allele: str,
                                  chrm_end_pos: int, chrm_path: pathlib.Path, one_hot_encoder, enformer) -> dict:
    """
    Given a DNA positive strand sequence of size Enformer_input_len, the function generates the mutated sequence and then
    runs Enformer twice: once on the input sequence and mutated sequence and once on the reversed-complement
    of the input sequence and mutated sequence.

    rel_sp: 0-based mutation start position in input sequence. Allowed values are 0 to len(ref_seq)-1.
    rel_sp: 0-based mutation end position in input sequence. Allowed values are 0 to len(ref_seq)-1.
    chrm_end_pos: 1-based chromosome position of the last NT in Enformer input.

    The output dictionary 'pred' has keys of False and True (corresponding to original sequence and reveresed
    complement sequence, respectively) and values that are dictionaries with keys 'ref' and 'mut', and
    values containning the (896, 5313) Enformer Human prediction of the reference and mutated sequences, respectively.
    That is:

        pref[rev][s] = Enformer (896,5313) human predicitons,
        
    where rev is either False or True, and s is either 'ref' or 'mut'.
    """
    assert len(ref_seq) == Enformer_input_len, f"Sequence length must be {Enformer_input_len:,}. It is {len(ref_seq)} !!"
    assert 0 <= rel_sp < len(ref_seq), f"{rel_sp=:,} is out of bound !!"
    assert 0 <= rel_ep < len(ref_seq), f"{rel_ep=:,} is out of bound !!"

    pred = {}

    # generate positive strand mutated sequence
    ref_seq_p = ref_seq.upper()
    mut_seq_p = generate_mut_variant(ref_seq_p, var_type, rel_sp, rel_ep, mut_allele).upper()
    
    # rev = False (i.e. positive sequence)
    # ------------------------------------
    c_ref_seq, c_mut_seq = ref_seq_p, mut_seq_p
    if len(c_mut_seq) != Enformer_input_len: # due to INS or DEL mutation
        c_mut_seq = adjust_seq_size(c_mut_seq, chrm_end_pos, str(chrm_path), Enformer_input_len, False)
    
    pred[False] = {
        'ref': np.squeeze(enformer_predict(enformer, one_hot_encoder(c_ref_seq))['human']),
        'mut': np.squeeze(enformer_predict(enformer, one_hot_encoder(c_mut_seq))['human'])
    }

    # rev = True  (i.e. reversed-complement sequence)
    # -----------------------------------------------
    c_ref_seq = reverse_complement(ref_seq_p)
    c_mut_seq = reverse_complement(mut_seq_p)
    if len(c_mut_seq) != Enformer_input_len: # due to INS or DEL mutation
        c_mut_seq = adjust_seq_size(c_mut_seq, chrm_end_pos, str(chrm_path), Enformer_input_len, True)

    pred[True] = {
        'ref': np.flipud(np.squeeze(enformer_predict(enformer, one_hot_encoder(c_ref_seq))['human'])),
        'mut': np.flipud(np.squeeze(enformer_predict(enformer, one_hot_encoder(c_mut_seq))['human']))
    }
    return pred


def is_in_Enformer_range(x: pd.Series, chrm_size: Optional[Dict] = None,
                         margin: int = 0,
                         start_pos_lbl: str = 'Start_Position',
                         tss_lbl: str = 'TSS',
                         chrm_lbl: str = 'Chromosome') -> bool:
    """
    Returns True [False] if the TSS is [is not] in the Enformer range (assuming mutation is in the middle of the range), 
    and [or] mutation is not [is] too close to the beginning or end of the chromosome.
    """
    s_pos = int(x[start_pos_lbl])
    flag = (np.abs(s_pos - int(x[tss_lbl])) <= (Enformer_predict_len/2)) and (s_pos > (Enformer_input_len//2 + margin))
    if chrm_size is not None:
        flag = flag and (np.abs(chrm_size[x[chrm_lbl]] - s_pos) > (Enformer_input_len//2 + margin))
    return flag


# metrics
def metric_wrapper_function(pred_f: dict, tss_bin: int, **kwargs: dict) -> dict:
    """TBD """
    #define parameters for the different metrics
    bin_options = ["mean", "meanTSS", "max", "maxTSS"]
    track_options = ["mean", "meanCAGE", "max", "maxCAGE"]
    strand_options = ["mean", "max"]

    # get the CAGE track indexes
    df_tracks = kwargs["df_tracks"]
    cage_indexes = df_tracks.query("track_type == 'CAGE'")['index'].tolist()
    
    # 0 -based 
    tss_bin = tss_bin - 1
  
    # calculate abs(relative difference) (common for all metrics) 
    diff_forward = abs((pred_f[False]["mut"] - pred_f[False]["ref"]) / (pred_f[False]["ref"])) 
    diff_reverse = abs((pred_f[True]["mut"] - pred_f[True]["ref"]) / (pred_f[True]["ref"]))
    
    #calculate a score for each metric:
    metrics_scores = {}
    
    for b_option in bin_options:
        for t_option in track_options:
            for s_option in strand_options:
                
                #merge bins
                temp_score_forward = merge_bins(diff_forward, tss_bin, b_option) 
                temp_score_reverse = merge_bins(diff_reverse, tss_bin, b_option)
                
                #merge tracks
                temp_score_forward = merge_tracks(temp_score_forward, cage_indexes, t_option) 
                temp_score_reverse = merge_tracks(temp_score_reverse, cage_indexes, t_option)
                
                #merge strands
                final_score = merge_strands(temp_score_forward, temp_score_reverse, s_option)
                
                #update the dictionary holding the results
                metric_name = f"b{b_option}_t{t_option}_s{s_option}"
                metrics_scores[metric_name] = final_score
                
    return(metrics_scores)


def merge_bins(data: np.ndarray, tss_bin: int, option: str) -> np.ndarray:
    """TBD"""
    start = tss_bin - 1 if tss_bin > 0 else 0
    if option == "mean":
        data = data.mean(axis = 0)
    elif option == "meanTSS":
        #data = data[tss_bin - 1: tss_bin + 2, :].mean(axis = 0)
        data = data[start : tss_bin + 2, :].mean(axis = 0)
    elif option == "max":
        data = data.max(axis = 0)
    elif option == "maxTSS":
        #data = data[tss_bin - 1: tss_bin + 2, :].max(axis = 0)
        data = data[start : tss_bin + 2, :].max(axis = 0)
        
    return(data)


def merge_tracks(data: np.ndarray, cage_indexes: list, option: str) -> float:
    """TBD"""
    if option == "mean":
        data = data.mean()
    elif option == "meanCAGE":
        data = data[cage_indexes].mean()
    elif option == "max":
        data = data.max()
    elif option == "maxCAGE":
        data = data[cage_indexes].max()
    else:
        return 0.0
    return(data)


def merge_strands(data_forward: float, data_reverse: float, option: str) -> float:
    """TBD"""
    if option == "mean":
        score = np.mean([data_forward, data_reverse])
    elif option == "max":
        score = np.max([data_forward, data_reverse])
    else:
        return 0.0
    return(score)


def run_enformer_on_dataframe(df: pd.DataFrame, store_path: pathlib.Path,
                              chrom_path: pathlib.Path, one_hot_encoder, enformer,
                              debug_interval: int = 50,
                              rslt_file_append_lbl: str = '',
                              metric_funcs: Optional[List[Callable[[dict, int, dict], np.float64]]] = None,
                              **metric_funcs_kwargs) -> None:
    """
    Runs Enformer on a given dataframe of mutations.

    The dataframe must include the following columns:
    'mut_id', 'Chromosome', 'Variant_Type', 'TSS', 'Transcript stable ID',
    'ref_allele', 'mut_allele', 'Start_Position', 'End_Position'

    Valid Variant_Type are SNP, DEL, and INS.

    The function saves the Enformer results for each mutation in a pickle file. The result is either
    the raw Enformer output or metrics of the raw Enformer output. The name
    of the saved pickle file is:
      mut_id.replace(':', '_').replace('-', 'X') + f"__{Transcript stable ID}_{TSS}.pkl

    The kwargs arguments are for the metric_funcs functions.
    """
    # the 0-based offset of the mutation start position in ref_seq_p below (mutation in the middle of the context)
    rel_sp = Enformer_input_len//2 - 1 + Enformer_input_len%2

    for i, (_, r) in enumerate(df.iterrows()):
        if (i > 0) and (i % debug_interval == 0):
            print(f"Processed {i} rows ...")

        mut_id, m_chrm, var_type, tss, tid, ref_allele, mut_allele, start_position, end_position = r[
            ['mut_id', 'Chromosome', 'Variant_Type', 'TSS', 'Transcript stable ID',
             'ref_allele', 'mut_allele', 'Start_Position', 'End_Position']
            ]
        start_position, end_position = int(start_position), int(end_position)

        #print(i, mut_id, var_type, tss, tid, ref_allele, mut_allele, start_position, end_position)

        # Enformer intervals
        m_chrm_file = chrom_path / get_chrm_file_name(m_chrm[3:])
        predt_start, predt_end, m_interval, in_m_start_pos, in_m_end_pos, ex_m_interval = enformer_intervals(start_position, m_chrm)

        # reference sequence
        ref_seq_p = pull_fasta_sequence(str(m_chrm_file), in_m_start_pos, in_m_end_pos).upper()
        if len(ref_seq_p) != Enformer_input_len:
            print(f"{i=}, {mut_id=} => reference sequence size ({len(ref_seq_p)}) different than required ({Enformer_input_len}). Maybe mutation located toward the begenning or end of the chromosome. Skipping ...")
            continue

        rel_ep = rel_sp + (end_position - start_position) # the 0-based offset of the mutation end position in ref_seq_p
        if var_type != 'INS': # in INS case, ref_allele=='-'
            if ref_seq_p[rel_sp : rel_sp+len(ref_allele)] != ref_allele:
                print(f"{i}. {mut_id=}, reference allele mismatch with chromosome position !!. Skipping ...")
                continue
        
        # run Enformer
        pred_f = enformer_predict_both_strands(ref_seq_p, var_type, rel_sp, rel_ep, mut_allele, in_m_end_pos, m_chrm_file, one_hot_encoder, enformer)
        
        # saving either the metrics or the raw output
        save_info = ('Enformer', pred_f) if metric_funcs is None else ('Enformer_metrics', {f.__name__: f(pred_f, pos2bin(predt_start, tss, rev=False)[0], metric_funcs_kwargs) for f in metric_funcs})

        # save to pickle
        
        # rslt_file = store_path / (mut_id.replace(':', '_').replace('-', 'X') +
        #                           (f"__{tid}_{tss}.pkl" if rslt_file_append_lbl == '' else f"__{tid}_{tss}_{rslt_file_append_lbl}.pkl")
        #                           )
        rslt_file = store_path / (mut_id.replace(':', '_').replace('-', 'X') + f"__{tid}_{tss}.pkl")
        with open(rslt_file, 'wb') as fp:
            pickle.dump(
                #{'mut_id': mut_id, 'TSS': tss, 'TID': tid, 'predt_start': predt_start, 'predt_end': predt_end, 'Enformer': pred_f}, fp, protocol=pickle.HIGHEST_PROTOCOL
                {'mut_id': mut_id, 'TSS': tss, 'TID': tid, 'predt_start': predt_start, 'predt_end': predt_end, save_info[0]: save_info[1]}, fp, protocol=pickle.HIGHEST_PROTOCOL
                )

        del pred_f  # dont think this is needed
        
def reverse_complement(s: str, complement: dict = None) -> str:
    """Performs reverse-complement of a sequence. Default is a DNA sequence."""
    if complement is None:
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    s_rev = s[::-1]
    lower = [b.islower() for b in list(s_rev)]
    bases = [complement.get(base, base) for base in list(s_rev.upper())]
    rev_compl = ''.join([b.lower() if l else b for l, b in zip(lower, bases)])
    return rev_compl

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
We want to use metrics that are specificaly relevant for the MDR1 gene. 
- Use only tracks that are relevant to tissues where MDR1 is expressed
- Since the TSS is far from the variants, we dont need to use it. 
we will create metrics that reflect our needs. 
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def metric_wrapper_function_mdr1(pred_f: dict, tss_bin, **kwargs: dict) -> dict:
    """TBD """
    #define parameters for the different metrics
    bin_options = ["mean"]
    track_options = ["mean", "meanCAGE", "meanMDR1", "meanCAGEMDR1"]
    strand_options = ["mean"]

    df_tracks = kwargs["df_tracks"]
    # get the CAGE track indexes
    cage_indexes = df_tracks.query("track_type == 'CAGE'")['index'].tolist()
    #get MDR1 track indexes
    mdr1_tissues = ["brain", "kidney", "colon", "liver"] #tissues where MDR1 is highly expressed
    mask = df_tracks['track_name'].str.contains('|'.join(mdr1_tissues), case=False) #"True" for tracks that are in any of these tissues
    mdr1_indexes = df_tracks.index[mask] 
    #get MDR1 & CAGE track indexes - intersect
    mdr1_cage_indexes = [ind for ind in cage_indexes if ind in mdr1_indexes]

    # calculate abs(relative difference) (common for all metrics) 
    diff_forward = abs((pred_f[False]["mut"] - pred_f[False]["ref"]) / (pred_f[False]["ref"])) 
    diff_reverse = abs((pred_f[True]["mut"] - pred_f[True]["ref"]) / (pred_f[True]["ref"]))
    
    #calculate a score for each metric:
    metrics_scores = {}
    
    for b_option in bin_options:
        for t_option in track_options:
            for s_option in strand_options:
                
                #merge bins
                temp_score_forward = merge_bins(diff_forward, tss_bin, b_option) 
                temp_score_reverse = merge_bins(diff_reverse, tss_bin, b_option)
                
                #merge tracks
                temp_score_forward = merge_tracks_mdr1(temp_score_forward, cage_indexes, mdr1_indexes, mdr1_cage_indexes, t_option) 
                temp_score_reverse = merge_tracks_mdr1(temp_score_reverse, cage_indexes, mdr1_indexes, mdr1_cage_indexes, t_option)
                
                #merge strands
                final_score = merge_strands(temp_score_forward, temp_score_reverse, s_option)
                
                #update the dictionary holding the results
                metric_name = f"b{b_option}_t{t_option}_s{s_option}"
                metrics_scores[metric_name] = final_score
                
    return(metrics_scores)

#TBD: continue to set the "metric_wrapper_function_mdr1" and the "merge_tracks_mdr1" to fit our case. 
#Then run on all variants, original and random. See if there is any metric where the results are significant. 

def merge_tracks_mdr1(data: np.ndarray, cage_indexes: list,  mdr1_indexes: list, mdr1_cage_indexes: list, option: str) -> float:
    """TBD"""
    if option == "mean":
        data = data.mean()
    elif option == "meanCAGE":
        data = data[cage_indexes].mean()
    elif option == "meanMDR1":
        data = data[mdr1_indexes].mean()
    elif option == "meanCAGEMDR1":
        data = data[mdr1_cage_indexes].mean()  
    else:
        return 0.0
    return(data)

