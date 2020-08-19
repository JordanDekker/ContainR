import re
from Bio.Seq import Seq
from app.core.exceptions.exceptions import NonDNAException

ALLOWED_EXTENSIONS = set(['fastq', 'csv', 'tsv', 'xml', 'out', 'zip', 'gz'])


def get_reverse_complement(input_seq):
    """Returns the reverse complement of a sequence.

    Args:
        input_seq: DNA or RNA sequence, must be a string.
    
    Returns:
        String: The reverse complement of input_seq.
    
    Raises:
        NonDNAException: If the sequence contains non-DNA characters.
        TypeError: If input_seq is not a string.

    """
    try:
        if check_dna(input_seq.upper()):
            return str(Seq(input_seq).reverse_complement())
        else:
            raise NonDNAException
    except TypeError:
        raise TypeError


def check_dna(seq, code=re.compile(r'[^ACGTN.]').search):
    """Checks if a sequence consists of DNA-characters.
    
    Args:
        seq: Sequence, must be a string.
        code: Regex pattern.
    
    Returns:
        True: If the sequence is pure DNA.
        False: If the sequence if not (pure) DNA.

    """
    return not bool(code(seq))


def allowed_file(filename):
    """Checks if a file is allowed by checking the extention.

    Args:
        filename: String of path to file.
    
    Returns:
        True: If the file is allowed.
        False: If the file is not allowed.
    
    """
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS
