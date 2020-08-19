import pandas as pd
import ast
import json
from flask import jsonify

def receive_json_ids(dataframe, jsondata, just_headers = False):
    """Filters a dataframe to create a json for visualisation.
    
    This function loads in a dataframe and open a json file. It takes
    the last entry from the json list and searches all
    the entries in a column this column is specified by the length of
    the json list. With these ids it creates a new
    dataframe with only the data of the specified ids this is returned
    as json so it can be used for visualisation.

    Args:
        dataframe: The main dataframe.
        jsondata: Data to filter on.
        just_headers: Boolean value. If True it will return lists of the headers.

    Returns:
        fw_headers: Lit of forward headers.
        rv_headers: Lit of reverse headers.
        response: Dictionary with data to display when a sunburst node is clicked.

    """

    dict_data = ast.literal_eval(jsondata)
    jsondict = {1: 'Kingdom', 2: 'Phylum', 3: 'Class', 4: 'Order', 5: 'Family', 6: 'Genus', 7: 'Species', 8 : 'Strain'}
    # this checks how long the jsondata is and from this it selects the correct Letter out of the jsondict #
    suffix = jsondict[len(dict_data)]

    # This selects the data which has the same name as the recieved jsondata
    fw_subset = dataframe[(dataframe["fw_"+ suffix] == dict_data[-1])] 
    rv_subset = dataframe[(dataframe["rv_"+suffix] == dict_data[-1])]

    # This is only used so that the columns can be easily renamed in something more generic so the append will merge the correct columns
    columns_rename = pd.DataFrame(columns=["bitscore", "identity", "length"])

    # Get the specified data
    fw_sideDf = fw_subset[["fw_bit", "fw_id", "fw_coverage_length"]]
    rv_sideDf = rv_subset[["rv_bit", "rv_id", "rv_coverage_length"]]

    # Get headers
    fw_headers = fw_subset.index.values.tolist()
    rv_headers = rv_subset.index.values.tolist()

    if just_headers:
        return fw_headers, rv_headers
        
    # Rename the columns
    fw_sideDf.columns = columns_rename.columns
    rv_sideDf.columns = columns_rename.columns
    # Combine the two dataframes in one since they have the same column names it will merge completly
    sideDf = fw_sideDf.append(rv_sideDf)
    # Count and group the different entries also convert them into a json
    count_id = sideDf.round(0).groupby(['identity']).size().to_json(orient='table')
    count_bit = sideDf.round(0).groupby(['bitscore']).size().to_json(orient='table')
    count_length = sideDf.round(0).groupby(['length']).size().to_json(orient='table')
    fw_seqs = fw_subset["fw_seq"].tolist()
    rv_seqs = rv_subset["rv_seq"].tolist()

    # Get taxonomy id's
    tax_ids = set([*fw_subset.fw_accession.tolist(), *rv_subset.rv_accession.tolist()])
    tax_len = len(tax_ids)
    if tax_len == 0:
        tax_id = "None"
    elif tax_len == 1:
        tax_id = list(tax_ids)[0]
    else:
        tax_id = "More"

    response = {
        "count_id":count_id,
        "count_bit": count_bit,
        "count_length": count_length,
        "node_name": dict_data[-1],
        "tax_id": str(tax_id),
        "fw_headers": fw_headers,
        "rv_headers": rv_headers,
        "fw_seqs": fw_seqs,
        "rv_seqs": rv_seqs
    }
    return jsonify(response)
