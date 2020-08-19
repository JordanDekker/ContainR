from os import path, sep

from flask import session
from app.core.containr.check_taxonomy import *
from app.core.containr.parse_blast_results import *
from app.core.containr.load_data import *
from app.core.containr.BLAST import *


def containr_main(ref_file, fastq_df, blastOptions):
    """Runs the blast script, taxonomic check script and updates the main dataframe.

    Args:
        ref_file: The reference file which contains UIDs of the 16srRNA genes of NCBI and the Taxonomic name linked to this UID.
        fastq_df: The main dataframe.
        blastOptions: Dictionary with options for the command line BLAST module.
    
    """
    session_id = session['id']
    main_df = fastq_df.get_dataframe()
    convert_to_fasta(main_df, 'data/' + session_id + '/fastain.fasta')
    make_blast_database()
    blast_me('data/' + session_id + '/fastain.fasta', 'data/' + session_id + '/blastout.xml', blastOptions)

    ref_dic = read_ref(ref_file)
    new_df = read_output(ref_dic, 'data/' + session_id + '/blastout.xml', main_df)
    fastq_df.set_df(new_df)
