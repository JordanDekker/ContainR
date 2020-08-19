import pandas as pd
from Bio import Align
from Bio.SubsMat import MatrixInfo
from sklearn.cluster import KMeans


def cluster_seq(main_df, prefix='', num_clusters=3):
    """Clusters the sequences to see similair sequences.

    Args: 
        main_df: The whole pandas dataframe-subset that you want to visualise.
        prefix: 'fw_' or 'rv_' optional prefix for if you want to cluster fw or rv sequences.
        num_clusters: The amount of clusters you want. Default value is 3.
    
    Returns: 
        main_df: A dataframe sorted by the clusters they're put in.

    """
    distance_df = pd.DataFrame(index=main_df.index, columns=main_df.index)

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = MatrixInfo.blosum62
    
    for main_df_row1 in main_df.index:
        for main_df_row2 in main_df.index[main_df_row1-1:]:
            sequence1 = main_df[prefix+'seq'][main_df_row1]
            sequence2 = main_df[prefix+'seq'][main_df_row2]
            score = aligner.score(sequence2, sequence1)
            distance_df.at[main_df_row1, main_df_row2] = score
            distance_df.at[main_df_row2, main_df_row1] = score
              
    km = KMeans(n_clusters=num_clusters, init='k-means++', n_init=10)
    km.fit(distance_df)
    clusters = km.fit_predict(distance_df)
    main_df[prefix+'cluster'] = clusters
    
    main_df = main_df.sort_values(by=[prefix+'cluster'])
    
    return main_df
