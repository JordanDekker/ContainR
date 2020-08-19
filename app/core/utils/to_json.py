import pandas as pd
import numpy as np
import collections
import json
import math
from statistics import median


def seq_length_to_json(df):
    """Convert sequence length distribution to JSON object.

    Args:
        df: DataFrame containing a subset of sequence length.
    
    Returns:
        String: A JSON-formatted string for visualization of sequence length distribution.

    """
    json_data = {}
    for column in list(df):
        json_values = []
        distribution = collections.Counter(df.loc[:, column].dropna().tolist())
        for key, value in distribution.items():
            json_values.append({"x": int(key), "y": value})
        json_data[column] = json_values
    return json.dumps(json_data)


def perc_count_to_json(df):
    """?"""
    df_count = pd.Series.to_frame(df)
    df_count.index.names = ["perc"]
    return df_count.to_json(orient="table")


def get_paired_percentage_to_json(df):
    """Get a JSON object containing the percentage of True and False paired reads.
    
    Args:
        df: Pandas DataFrame.
    
    Returns:
        String: A JSON-formatted string.

    """
    df = df.loc[df["paired_flag"] == False]
    df = df if len(df.index) > 0 else df.loc[df["paired_flag"] == True]
    json_paired_data = []
    paired_seqs = df.groupby(["paired"]).size()
    if paired_seqs.count() < 2:
        json_paired_data.append({"name": "True", "y": 100})
        json_paired_data.append({"name": "False", "y": 0})
        return json.dumps(json_paired_data)
    else:
        paired_seq_number = paired_seqs.get_values()
        true_values = paired_seq_number[1]
        false_values = paired_seq_number[0]
        total = true_values + false_values
        true_values_percentage = round((true_values/total)*100, 3)
        false_values_percentage = round((false_values/total)*100, 3)
        json_paired_data.append({"name": "True", "y": true_values_percentage})
        json_paired_data.append({"name": "False", "y": false_values_percentage})
        return json.dumps(json_paired_data)


def nucleotide_percentages_to_json(df):
    """Calculate box plot values from nucleotide percentages.
    
    Args:
        df: A pandas dataframe containing a subset of nucleotide percentages.
    
    Returns:
        json_data: A JSON-formatted list for visualization of nucleotide percentages.

    """
    json_data = []

    for pair in enumerate(df.iteritems()):
        data = list(pair[1][1].values)
        data = sorted(data)

        # Calculate quartiles
        q1, q3 = np.percentile(data,[25,75])
        q2 = median(data)
        
        # Calculate IQR (Interquartile Range)
        iqr = q3 - q1
        
        # Define bounds
        lower_bound = q1 -(1.5 * iqr) 
        upper_bound = q3 +(1.5 * iqr) 
        
        # Get outliers
        outliers = [i for i in data if i < lower_bound or i > upper_bound]
        non_outliers = [i for i in data if i >= lower_bound and i <= upper_bound]

        # Calculate total reads in every quarter. This includes outliers.
        # TODO Think about where to use >= or > and < or <=.
        quarter1 = len([i for i in data if i < q1])
        quarter2 = len([i for i in data if i >= q1 and i < q2])
        quarter3 = len([i for i in data if i >= q2 and i < q3])
        quarter4 = len([i for i in data if i >= q3])
        
        # Min and max within bounds, median, Q's.
        # The order of these matters for the visualisation in the graph.
        box_values = [min(non_outliers), q1, q3, max(non_outliers), q2]
        json_data.append({'x':pair[0],'label':pair[1][0],'y':box_values, 'outliers': outliers, 'quarter1': quarter1, 'quarter2': quarter2, 'quarter3': quarter3, 'quarter4': quarter4})

    return json_data
