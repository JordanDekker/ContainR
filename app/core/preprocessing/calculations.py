import numpy as np
import pandas as pd

import multiprocessing as mp
from Bio import Align


def get_paired_reads(df):
    """Check if the forward and reverse complement are both present.
    
    Adds a column to the dataframe called "paired" which has the value True if both sequences are present.
    Args:
        df: Dataframe which has the columns "fw_seq" and "rvc_seq".
    
    Returns:
        Dataframe with the column "paired" which has True if both sequences are present.

    """
    df['paired'] = df[["fw_seq", "rv_seq"]].notna().all(axis=1)


def get_nucleotide_percentages(df, prefix):
    """Calculate the percentages of nucleotides in a sequence.
    
    Adds columns to the dataframe with percentages of nucleotides.

    Args:
        df: The main df.
        prefix: Prefix stating the strand, either 'fw_' or 'rvc_', must be a string.
    
    Raises:
        NonDNAException: If the sequence contains non-DNA characters.
        AttributeError: If input_seq has no attribute 'upper' because it is not a string.
        TypeError: If input_seq or prefix is not a string.
        ZeroDivisionError: If input_seq is an empty string.

    """
    try:
        df[prefix + 'A_perc'] = df[prefix + 'seq'].apply(lambda row: round(row.upper().count('A') / len(row) * 100, 4) if pd.notnull(row) else np.nan)
        df[prefix + 'T_perc'] = df[prefix + 'seq'].apply(lambda row: round(row.upper().count('T') / len(row) * 100, 4) if pd.notnull(row) else np.nan)
        df[prefix + 'G_perc'] = df[prefix + 'seq'].apply(lambda row: round(row.upper().count('G') / len(row) * 100, 4) if pd.notnull(row) else np.nan)
        df[prefix + 'C_perc'] = df[prefix + 'seq'].apply(lambda row: round(row.upper().count('C') / len(row) * 100, 4) if pd.notnull(row) else np.nan)
    except TypeError:
        raise TypeError
    except KeyError:
        raise KeyError


def get_sequence_length(df, doReturn = False):
    """Calculates length for all sequences and adds this to the column x_seq_length.
    
    Adds columns to the dataframe with the total length of sequences.

    Args:
        df: The dataframe which to insert the values into.
    
    Raises:
        TypeError: raised if the sequences aren't strings.
        KeyError: raised if the label does not exist in the dataframe.

    """
    try:
        df['fw_seq_length'] = pd.array(df['fw_seq'].apply(lambda fw_len: len(fw_len) if pd.notnull(fw_len) else np.nan), dtype="Int64")
        df['rv_seq_length'] = pd.array(df['rv_seq'].apply(lambda rv_len: len(rv_len) if pd.notnull(rv_len) else np.nan), dtype="Int64")
        if doReturn:
            return df
    except TypeError:
        raise TypeError
    except KeyError:
        raise KeyError

def make_alignment(df,matchScore,mismatchScore, open_gap_score, extend_gap_score):
    tempdf = df.loc[:,["fw_seq","rvc_seq"]]
    tempdf["matchScore"] = matchScore
    tempdf["gapScore"] = open_gap_score
    tempdf["extend_gap_score"] = extend_gap_score
    tempdf["mismatchScore"] = mismatchScore
    alignment = list(map(calculate_identity, tempdf.values.tolist()))
    return(add_alignment(df, alignment,matchScore))

def calculate_identity(pairList):
    try:
        if pd.isna(pairList[0]) or pd.isna(pairList[1]):
            raise IndexError
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = pairList[2]
        aligner.open_gap_score = pairList[3]
        aligner.extend_gap_score = pairList[4]
        aligner.mismatch_score = pairList[5]
        alignment = aligner.align(pairList[0],pairList[1])
        return str(alignment.score)+";"+str(alignment[0]).replace("\n",";")
    except IndexError:
        return "0;;;"
  
def add_alignment(df,alingment,matchScore):
    indexDict = {0 : "overlap_identity_perc",1 : "fw_alignment", 2 : "aligmnent_intersection", 3 : "rvc_alignment"}
    for key in indexDict.keys():
        df[indexDict[key]] = list(map(lambda i: i.split(";")[key] if i is not None else i,alingment))
    df["overlap_identity_perc"] = pd.to_numeric(df["overlap_identity_perc"])
    df["max_possible_score"] = matchScore
    df["max_possible_score"] = list(map(get_max_score, df.loc[:,["fw_seq_length","rv_seq_length","max_possible_score"]].values.tolist()))
    df["overlap_identity_perc"] = list(map(lambda x: x[0]/x[1]*100 if x is not None else None, df.loc[:,["overlap_identity_perc","max_possible_score"]].values.tolist()))
    df = df.drop(columns = ["max_possible_score"])
    return df

def get_max_score(row):
    if None not in row:
        return min([row[0]*row[2],row[1]*row[2]]) 
    else:
        return None
