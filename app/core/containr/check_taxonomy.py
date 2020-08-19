from Bio import Entrez
from collections import defaultdict
import time
from flask import session


def read_file(ref_file):
    """Reads a reference file and converts it to a dictionary.
    
    Args:
        ref_file: String of the path to the file that contains the reference data.

    Returns: 
        gene_UID_dict: Taxonomic_UID_dic, this is an dictionary where the key is the UID and the value the taxonomic information.

    """
    gene_UID_dict = {'':[]}
    with open(ref_file,'r') as file:
        for line in file:
            ref_record = line.replace("\n", "").split(",")
            gene_UID_dict[ref_record[0]] = ref_record
        # deletes records with empty values #
        del gene_UID_dict['']
    return gene_UID_dict
    

def get_taxonomic_ID(gene_UID_dict):
    """Uses the gene UID to get the taxonomic UID.

    Args: 
        gene_UID_dict: A dictionary where the gene_UID are keys and the taxonomic name from the file is the value.
    
    Returns: taxonomic_UID_dic: A dictionary where the taxonomy UIDs are keys and the taxonomic name is the value.

    """
    Entrez.email = "thijs.schoppema@gmail.com"
    gene_UID_list = list(gene_UID_dict.keys())
    taxonomic_UID_dict = defaultdict(list)
    while len(gene_UID_list) != 0:
        time.sleep(5)
        NCBI_data = Entrez.esummary(db="nucleotide", id=','.join(gene_UID_list[0:4000])) # ncbi max cap
        loaded_NCBI_data = Entrez.read(NCBI_data)
        for line in loaded_NCBI_data:
            taxonomic_UID_dict[str(line['TaxId'])].append(gene_UID_dict[gene_UID_list.pop(0)])
    return taxonomic_UID_dict


def get_taxonomy(taxonomic_UID_dict):
    """Uses the gained taxonomy UIDs to querry the taxonomy db and get the full taxonomic path.

    Args: 
        taxonomic_UID_dict: A dictionary where the taxonomy UIDs are keys and the taxonomic name is the value.

    Returns:
        taxonomic_string: A string of the results with the fixed taxonomy.

    """
    Entrez.email = "thijs.schoppema@gmail.com"
    id_list = list(taxonomic_UID_dict.keys())
    taxonomic_string = ''
    taxonomic_set = set(['superkingdom', 'phylum','class','order','family','genus','species'])
    while len(id_list) != 0:
        time.sleep(5)

        NCBI_data = Entrez.efetch(db="taxonomy", id=','.join(id_list[0:4000]))     #arround the max input NCBI allows
        loaded_NCBI_data = Entrez.read(NCBI_data)
        for record in loaded_NCBI_data:
            tax_record_dict = {'superkingdom':'u-k', 'phylum':'u-p', 'class':'u-c', 'order':'u-o', 'family':'u-f', 'genus':'u-g', 'species':'u-s'}
            for taxonomic_data in record['LineageEx']:
                tax_record_dict[taxonomic_data['Rank']] = taxonomic_data['ScientificName']
            tax_record_dict['species'] = tax_record_dict['genus'][0] + '_' + record['ScientificName'].split(' ')[1]    #We only want the species from this
            record_taxonomic_set = set(tax_record_dict.keys())
            record_taxonomic_set = record_taxonomic_set.difference(taxonomic_set)
            for taxonomic_data in record_taxonomic_set:
                del tax_record_dict[taxonomic_data]
            ref_lijst = taxonomic_UID_dict[id_list[0]]       #the first id of the list/the one used right now

            for i, ref_record in enumerate(ref_lijst):
                ref_record[2] = '; '.join(list(tax_record_dict.values()))  #ref record[2] is the old taxonomic value thats been overwritten
                taxonomic_string = taxonomic_string + ','.join(ref_record) + '\n'
            del id_list[0]  #remove the id from the list
    return taxonomic_string


def save_result(new_ref_file, ref_file):
    """Overwrites the reference set with the corrected taxonomy.
    
    Args:
        new_ref_file: String of path to new reference file.
        ref_file: String of path to old reference file.
        
    """
    with open(ref_file,'w+') as file:
        file.write(new_ref_file)
