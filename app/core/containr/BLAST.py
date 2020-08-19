from multiprocessing import cpu_count
from subprocess import call

from flask import session


def make_blast_database():
    """Creates a BLAST database.

    Creates database using the command line on which the blast+ module can BLAST more efficiently.
    It requires a fasta file with the data over which you want to blast.

    """
    session_id = session['id']
    cmd = "makeblastdb -in app/core/containr/Ref_files/ref_file.fasta -dbtype nucl -out data/" + session_id + \
          "/blastndatabase"" -parse_seqids"
    call(cmd, shell=True)


def blast_me(fasta_file, blast_file, blastOptions):
    """Executes a BLAST.

    Use a fasta input file to blast over the blastn database. 
    The output is an XML file named blast_out.xml
    
    Args:
        fasta_file: Path to the fasta file.
        blast_file: Path to save BLAST output to.
        blastOptions: Dictionary with options for the command line BLAST module.

    """
    session_id = session['id']
    if blastOptions["use_options"]:
        cmd = f"blastn -query {fasta_file} -db data/{session_id}/blastndatabase -out {blast_file}" + \
            f" -outfmt 5 -max_target_seqs 1 -num_threads {str(cpu_count())} -perc_identity {blastOptions['minIdent']}" + \
            f" -reward {blastOptions['reward']} -penalty {blastOptions['penalty']} -gapopen {blastOptions['gapOpen']} -gapextend {blastOptions['gapExtend']} -word_size 50"
    else:
        cmd = f"blastn -query {fasta_file} -db data/{session_id}/blastndatabase -out {blast_file}" + \
        f" -outfmt 5 -max_target_seqs 1 -num_threads {str(cpu_count())} -word_size 50"

    result = call(cmd, shell=True)
