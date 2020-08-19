import pandas as pd

DEFAULT_INIT_FOLDERS = ['OUTPUT', 'ANNOTATION']
METADATA_FILE_PATH = 'metadata.txt'


def load_input(input_path: str):
    """Converts TSV file to pandas dataframe.
    
    Loads input data from a tab-separated values file and converts it to a pandas dataframe.
    Excludes non-flagged items and returns only columns containing sequences.

    Args:
        input_path: String of path to the input file.
    
    Returns:
        df: Dataframe with forward and reverse complement sequences.
    
    Raises:
        KeyError: If requested key (e.g. sequence) is absent and can't be loaded.
        FileNotFoundError: If file_path doesn't refer to an existing file.
        ValueError: If an incorrect object type is used.

    """
    try:
        df = pd.read_table(input_path, sep='\t', header='infer', index_col=0, comment='#')
        return df[~df.flagged]
    except KeyError:
        raise KeyError
    except FileNotFoundError:
        raise FileNotFoundError
    except ValueError:
        raise ValueError


def convert_to_fasta(df: pd.DataFrame, output_path: str):
    """Convert sequences and headers in DataFrame to FASTA-format.
    
    Include postfixes '/1' for forward- and '/2' for reverse complement sequences in FASTA header.
    Saves output to  file with sequences in FASTA-format.

    Args:
        df: Dataframe containing sequences (columns) and headers (row indices), type must be pd.DataFrame.
        output_path: Path to output file, type must be str.

    Raises:
        KeyError: If requested key (e.g. sequence) is absent and can't be loaded.
        ValueError: If an incorrect object type is used.

    """
    try:
        content = ''
        with open(output_path, 'w') as f:
            filtered_df = df[~df.flagged]
            for index, row in filtered_df.iterrows():
                if isinstance(row['fw_seq'], str):
                    content += '>' + index + ' /1\n' + row['fw_seq'] + '\n'
                if isinstance(row['rvc_seq'], str):
                    content += '>' + index + ' /2\n' + row['rvc_seq'] + '\n'
            f.write(content)
    except KeyError:
        raise KeyError
    except ValueError:
        raise ValueError
