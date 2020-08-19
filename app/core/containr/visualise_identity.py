import math


def perform_cutoff(df, steps, prefix, col = 'id'):
    """

    Args:
        df: The main dataframe that has to be changed/filtered
        steps: The steps in which you want to filter: 1% means that each 1% 1 taxonomic rang will be removed. (???)
        prefix: 'rv_' or 'fw_' to select either the forward or the reverse reads.
        col: The name of the column in which you want to filter. Length and bit score aren't supported yet.

    Returns:
        df: The filtered dataframe.
        
    """
    tax_list = [prefix+'Kingdom', prefix+'Phylum', prefix+'Class', prefix+'Order', prefix+'Family', prefix+'Genus', prefix+'Species', prefix + 'Strain']
    for index, row in df.iterrows():
        if not math.isnan(row[prefix+col]) and row[prefix+'full_tax'] != None:
            cutoff = math.floor((100-row[prefix+col])/steps)
            df.at[index, prefix+'tax_cut_off'] = cutoff
            taxa = row[prefix + 'full_tax'].split('; ')
            cutoff = cutoff - (7 - len(taxa))
            if cutoff > 0:
                if cutoff >= 7:
                    df.at[index, prefix+'full_tax'] = None
                else:
                    taxa = taxa[0:-cutoff]
                    df.at[index, prefix+'full_tax'] = '; '.join(taxa)
                if cutoff > 7:
                    cutoff = 7
                for taxa_rang in tax_list[7-cutoff:]:
                    df.at[index, taxa_rang] = None
        
    return df
