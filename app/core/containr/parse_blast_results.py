import pandas as pd
from Bio.Blast import NCBIXML
from flask import session

def read_output(ref_dic, blast_xml_file, main_df):
    """Uses the Bio.Python parser to parse the output file of the BLAST.
    
    Args:
        ref_dic: This is an dictionary where the key is the UID and the value the taxonomy.
        result_file: This is the path to the file that contains the blast results in XML format.
        df: The df that is based upon the TSV input file.
    
    Returns:
        main_df: The original df with the blast results added in its own columns.
    
    """
    try:
        with open(blast_xml_file, 'r') as blast_xml:

            lijst = [{}, {}]
            # NCBIXML.parse makes a iterater within are objects with the data from the blast #
            blast_records = NCBIXML.parse(blast_xml)

            for record in blast_records:
                
                blast_record_alignment = record.alignments
                header = record.query

                if len(blast_record_alignment) != 0:
                    # selecting the wanted data from the record #
                    bit_score = blast_record_alignment[0].hsps[0].bits
                    length = blast_record_alignment[0].hsps[0].align_length
                    identity = int(100 / length * blast_record_alignment[0].hsps[0].identities)
                    accession = blast_record_alignment[0].hit_id
                    if '|' in accession:
                        accession = accession.split('|')[1]
                    taxonomic_name = ref_dic[accession]

                    if header[0] != '@':
                        header = '@'+header

                    # splitting the data in Forward sequences and Reverse sequences #

                    if '/1' in header:
                        lijst[0][header.replace(' /1', '')] = (
                            [accession, taxonomic_name, identity, bit_score, length] + taxonomic_name.split('; '))
                    else:
                        lijst[1][header.replace(' /2', '')] = (
                            [accession, taxonomic_name, identity, bit_score, length] + taxonomic_name.split('; '))

            for i, prefix in enumerate(['fw_', 'rv_']):
                #overwrite values even if indices dont match
                
                if prefix+'full_tax' in main_df.columns:
                    for column_name in ['accession', 'full_tax','id','bit','coverage_length','Kingdom','Phylum','Class','Order','Family','Genus','Species','Strain']:
                        del main_df[prefix+column_name]

                # transfer the list(lijst) with data to a dataframe #        
                taxonomic_df = pd.DataFrame.from_dict(lijst[i],
                                                columns=[prefix + 'accession', prefix + 'full_tax', prefix + 'id',
                                                         prefix + 'bit', prefix + 'coverage_length',
                                                         prefix + 'Kingdom', prefix + 'Phylum', prefix + 'Class', prefix + 'Order',
                                                         prefix + 'Family', prefix + 'Genus', prefix + 'Species', prefix + 'Strain'], orient='index')
                # combine the new dataframe with the main dataframe and check if the main_df isn't empty else update it
                
                    
                if len(main_df) == 0:
                    main_df.update(taxonomic_df) 
                
                temp_main_df = pd.merge(main_df,taxonomic_df, left_index=True, right_index=True, how="outer")

                dtypes = main_df.dtypes.combine_first(taxonomic_df.dtypes)
                # Make sure the columns are the correct dtype
                for key, value in dtypes.iteritems():
                    try:
                        temp_main_df[key] = temp_main_df[key].astype(value)
                    except ValueError:
                        if value in [int,float]:
                            temp_main_df[key] = pd.to_numeric(temp_main_df[key], errors = "coerce")
                        else:
                            raise ValueError

                    
                main_df = temp_main_df
                for column in [prefix + 'Kingdom', prefix + 'Phylum', prefix + 'Class', prefix + 'Order',
                        prefix + 'Family', prefix + 'Genus', prefix + 'Species', prefix + 'Strain']:
                        main_df[column] = list(map(fill_unknown_single, main_df[[column, prefix+"flagged"]].values.tolist()))
                del temp_main_df

        for strand in ["fw","rv"]:
            main_df[strand+"_full_tax"] = list(map(fill_unknown, main_df[[strand+"_full_tax",strand+"_flagged"]].values.tolist()))
            main_df[strand+"_Species"] = list(map(lambda x: x[0].split("; ")[-1] if len(str(x[0]).split("; ")) == 7 else x[1], main_df[[strand+"_full_tax",strand+"_Species"]].values.tolist()))
        return main_df
    except FileNotFoundError:
        raise FileNotFoundError


def fill_unknown(info):
    """Fill in the empty taxonomy spots in the full_tax column.

    Args:
        info: List containing taxonomy data.
    
    Returns:
        String with taxonomy data separated by semicolons.

    """
    tax = info[0]
    if info[1]: # info[1] is the flag column
        return tax
    if tax is None or tax == "" or pd.isnull(tax):
        return "unknown"
    else:
        tax = ["unknown" if x.startswith("u-") or x.startswith("unclassified") else x for x in tax.split("; ")]
        if len(tax) == 7:
            if "_" in list(tax[6]):
                tax[6] = tax[5]+" "+tax[6].split("_")[1]
        return "; ".join(tax)
        

def fill_unknown_single(info):
    """Fills in the empty tax spots per column.
    
    Args:
        info: ?

    Returns:
        item: ?
    
    """
    item = info[0]
    if info[1]: # check if flagged
        return item

    if item == "" or pd.isna(item) or pd.isnull(item):
        return "unknown"
    if item.startswith("u-") or item.startswith("unclassified"):
        return "unknown"
    else:
        return item


def read_ref(ref_file):
    """Reads the reference file and converts it to a dictionary.

    Args:
        ref_file: This is the path to the file that contains the reference data.
    Returns:
        ref_dic: This is an dictionary where the key is the UID and the value the taxonomy.

    """
    try:
        ref_dic = {}
        with open(ref_file, 'r') as file:
            for line in file:
                # Results from split(): index 0 is the first value that we take (UID), the 3 is the last indice,
                # and the 2 is the steps of what location we pick, so when use this we get the 0 and 2 indices that
                # stand for the UID and taxonomy.
                accession_code, taxonomic_name = line.replace('\n', '').split(',')[0:3:2]
                ref_dic[accession_code] = taxonomic_name
        return ref_dic
    except FileNotFoundError:
        raise FileNotFoundError


def convert_results(main_df):
    """Takes the dataframe and converts it into an json that works with the d3 sunburst visualisation.
    
    Args:
        main_df: Main dataframe object.
    
    Returns:
        ref_dic: A dictionary where the key is the UID and the value the taxonomy.

    """
    # Only select the full taxonomy and group them by name.
    taxonomy_df = pd.concat([main_df[main_df["fw_visual_flag"] == False]['fw_full_tax'], main_df[main_df["rv_visual_flag"] == False]['rv_full_tax']], ignore_index = True)
    # Drop the empty values.
    # TODO: These values are not usueless so find a way to use them.
    taxonomy_df = taxonomy_df.dropna()
    
    # Count the repeating stuff and sort them to make sure that double values won't get thrown.
    taxonomy_count = taxonomy_df.value_counts()
    taxonomy_count.sort_index(inplace=True)
    
    # Create the root.
    root = {'name':'flare', 'children':[]}
    # Loop over the species.
    for idx, species in enumerate(taxonomy_count):
        # Determine the sequence and size after wards split the sequence/taxonomy and reset the current_noce to the root/base.
        taxonomy_sequence = taxonomy_count.index[idx]
        size = int(species)        
        parts = taxonomy_sequence.split("; ") 
        current_node = root
        
        # Loop over the parts of the taxonomy and meanwhile traverse thew the dictionary/json accordingly.
        for i, node_name in enumerate(parts):
            children = current_node["children"]
            
            # Check if youre not at the end of the taxonomy chain.
            if i < (len(parts)-1):
                # Param to make sure that there will be no double childs/taxa in the sunburst.
                missing_child = True
                
                # Check the child nodes of the current taxa if the child is found make a record of it, set the missing_child to false because you found it.
                # After ward break out of the loop because you really dont want to loop over 200 extra values.
                for child in children:
                    if child['name'] == node_name:
                        child_node = child
                        missing_child = False
                        break
                
                # Check if we've found our missing taxa, if not add it to the dictionary.
                if missing_child:
                    child_node = {'name':node_name, 'children':[]}
                    children.append(child_node)
                
                # Set the current node in the dictionary as the current taxa.
                current_node = child_node
            
            # If we are at the end of the taxa make a specific one with size, name and childeren for options and add it to the dictionary.
            else:
                child_node = {'name':node_name, 'size':size, 'children':[]}
                children.append(child_node)
    return root


def main_parse_results(ref_file, blast_output, main_df):
    """Checks files and calls other functions.
    
    Args:
        ref_file: String of the path of the reference file.
        blast_output: The data recieved from the blastn funtion.
        main_df: The main dataframe.
    
    Return:
        main_df: The main dataframe.
        
    """
    session_id = session['id']
    try:
        ref_dic = read_ref(ref_file)
        main_df = read_output(ref_dic, blast_output, main_df)
        return main_df
    except FileNotFoundError:
        with open('data/' + session_id + "errors.txt", 'w+') as file:
            file.write("one of the files arn't pressent in the input folder")
    except RecursionError:
        with open('data/' + session_id + "errors.txt", 'w+') as file:
            file.write("the reference set was invalid")
