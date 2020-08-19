import re
import os
import json
import pandas as pd
import zipfile
from app.core.objects.FastQData import FastQData
from app.core.objects.FastQDataframe import FastQDataframe
from app.core.exceptions.exceptions import MissingDataException
from app.core.utils.preprocess_utils import get_reverse_complement
from app.core.preprocessing.calculations import get_nucleotide_percentages, get_paired_reads, get_sequence_length

PATTERN = re.compile("(^@[a-zA-Z0-9_].*\n)(.*\n)(\+\n)(.+($\n|.))", flags=re.MULTILINE)


def preprocess_fastq_files(og_fw_file, og_rv_file, uuid, options = None, flags = None):
    """Parse (zipped) FASTQ files into dataframe, extend this dataframe and save it to a Python pickle.

    Args:
        og_fw_file: Original forward FASTQ file.
        og_rv_file: Original reverse FASTQ file.
        Unique identifier, used as directory when saving the dataframe.

    """
    fw_file, rv_file = handle_zip(og_fw_file, og_rv_file, uuid)
    fastq_df = initialize_dataframe(uuid, fw_file, rv_file)
    extend_dataframe(fastq_df.get_dataframe())
    fastq_df.flag_any()
    fastq_df.set_preprocessOptions()
    if options is not None:
        configureOptions(options, fastq_df)
        if flags is not None:
            setFlags(flags, fastq_df)
    fastq_df.get_dataframe().to_csv("test.tsv",sep = "\t")
    fastq_df.to_pickle(path='data/'+uuid+'/', filename='pickle')


def preprocess_tsv_file(tsv_file, uuid):
    """"Parse tsv file into dataframe, extends the fastqdataframe object and returns the chosen options from the previous session.

    Args:
        tsv_file: Path to the TSV file.
        uuid: Unique identifier, used as directory when saving the dataframe.

    Returns:
        False: False boolean value.

    """
    options, data = splitTSV(tsv_file)
    with open(tsv_file,"w") as file:
        file.write(data)
    fw_data, rv_data, meta_data, flags = extractDataframe(tsv_file)
    fw_file = makeFastQ(fw_data,"fw",uuid)
    rv_file = makeFastQ(rv_data,"rv",uuid)
    preprocess_fastq_files(fw_file, rv_file, uuid, options, flags)

    if len(meta_data.columns) != 0:
        return addMetaData(uuid, meta_data)

    return False


def initialize_dataframe(uuid, fw_file, rv_file):
    """Parse the FASTQ files into FastQData objects and create a FastQDataframe using these FastQData objects.
    
    Args:
        uuid: Unique identifier.
        fw_file: Forward FASTQ file.
        rv_file: Reverse FASTQ file.

    Returns:
        fastaq_df: The initialized FastQDataframe.

    """
    fw_data = parse_fastq(os.path.abspath(fw_file), False)
    rvc_data = parse_fastq(os.path.abspath(rv_file), True)

    fastq_df = FastQDataframe(fw_data=fw_data, rvc_data=rvc_data, df_id=uuid)
    return fastq_df


def initialize_empty_dataframe(uuid):
    """Create and save an empty dataframe.

    Args:
        uuid: Unique identifier.

    """
    fastq_df = FastQDataframe(df_id=uuid, from_blast=True)
    fastq_df.to_pickle(path='data/'+uuid+'/', filename='pickle')


def extend_dataframe(fastq_df):
    """Extend the dataframe.

    Args:
        fastq_df: FastQDataframe to extend

    """
    get_paired_reads(fastq_df)
    get_nucleotide_percentages(fastq_df, "fw_")
    get_nucleotide_percentages(fastq_df, "rv_")
    get_sequence_length(fastq_df)


def parse_fastq(filename, is_rev):
    """Parse a FASTQ file into a FastQData object.

    Args:
        filename: Name of the file.
        is_rev: boolean to check if current file is a reverse FASTQ file.
    
    Returns:
        FASTQ files parsed into a FastQData object.
    """
    index_list = []
    index_info_list = []
    sequence_list = []
    complement_list = []
    q_score_list = []

    with open(os.path.abspath(filename)) as file_object:
        matches = re.findall(PATTERN, file_object.read())

        for match in matches:
            try:
                dna_sequence = match[1].strip()
                new_dna_sequence = re.sub(r'[^.ACGTN]', "N", dna_sequence)

                q_score_list.append(match[3].strip())
                sequence_list.append(new_dna_sequence)

                if is_rev:
                    complement_list.append(get_reverse_complement(new_dna_sequence))
                if " " in match[0]:
                    index_list.append(match[0].strip().split()[0])
                    index_info_list.append(match[0].strip().split()[1])
                else:
                    index_list.append(match[0].strip()[:-2])
            except Exception:
                raise MissingDataException("FASTQ file is incomplete")

    return FastQData(index_list, index_info_list, sequence_list, complement_list, q_score_list, is_rev)


def handle_zip(fw_file, rv_file, uuid):
    """Unzip forward and reverse FASTQ files if zipped.
    
    Args:
        fw_file: Original (zipped) forward FASTQ file.
        rv_file: Original (zipped) reverse FASTQ file.
        uuid: Unique identifier.
    
    Returns:
        fw_file: Unzipped forward files.
        rv_file: Unzipped reverse files.

    """
    if fw_file.rsplit('.', 1)[1] in ['zip', 'gz']:
        with zipfile.ZipFile(fw_file, 'r') as fw_zip:
            fw_zip.extractall("data/"+uuid+"/fw_file/")
            for root, dirs, files in os.walk("data/"+uuid+"/fw_file/"):
                for file in files:
                    fw_file = os.path.join(root, file)

    if rv_file.rsplit('.', 1)[1] in ['zip', 'gz']:
        with zipfile.ZipFile(rv_file, 'r') as fw_zip:
            fw_zip.extractall("data/"+uuid+"/rv_file/")
            for root, dirs, files in os.walk("data/"+uuid+"/rv_file/"):
                for file in files:
                    rv_file = os.path.join(root, file)
    return fw_file, rv_file


def extractDataframe(tsv_file):
    """Parse the tsv to a pandas dataframe and returns the data necessary to create fastQ files.

    Args:
        tsv_file: String of path to TSV file.

    Returns:
        List: A list of dataframes.

    """
    df = pd.read_csv(tsv_file,sep = "\t", header=0,index_col=0)
    return (df[['original_fw_seq','original_fw_seq_score']].fillna(value=""),
           df[['original_rv_seq','original_rv_seq_score']].fillna(value=""),
           df.drop(['original_fw_seq', 'original_fw_seq_score', 'fw_seq',
                    'fw_seq_score', 'fw_seq_length', 'fw_A_perc', 'fw_T_perc',
                    'fw_G_perc', 'fw_C_perc', 'fw_Q-score', 'original_rv_seq',
                    'original_rvc_seq', 'original_rv_seq_score', 'rv_seq', 'rvc_seq',
                    'rv_seq_score', 'rv_seq_length', 'rv_A_perc', 'rv_T_perc', 'rv_G_perc',
                    'rv_C_perc', 'rv_Q-score', 'flagged', 'fw_flagged', 'rv_flagged', 'paired_flag', 'fw_a_perc_flag',
                    'fw_t_perc_flag', 'fw_g_perc_flag', 'fw_c_perc_flag', 'rv_a_perc_flag',
                    'rv_t_perc_flag', 'rv_g_perc_flag', 'rv_c_perc_flag', 'fw_seq_len_flag', "original_contains_unknown_symbols",
                    'rv_seq_len_flag', 'identity_flag', 'fw_quality_flag', 'rv_quality_flag', 'paired'], axis = 1),
            df[['flagged', 'fw_flagged', 'rv_flagged', 'paired_flag', 'fw_a_perc_flag',
                    'fw_t_perc_flag', 'fw_g_perc_flag', 'fw_c_perc_flag', 'rv_a_perc_flag',
                    'rv_t_perc_flag', 'rv_g_perc_flag', 'rv_c_perc_flag', 'fw_seq_len_flag', "original_contains_unknown_symbols",
                    'rv_seq_len_flag', 'identity_flag', 'fw_quality_flag', 'rv_quality_flag', 'paired']])


def makeFastQ(fastqDF, strand, sessionID):
    """Create fastQ file out of a dataframe.

    Args:
        fastqDF: Pandas dataframe with sequences and scores.
        strand: Id if the data is about forward or reverse reads.
        sessionID: Unique identifier, used as directory when saving the fastq file.

    Returns:
        path: String of path to generated fastQ file.

    """
    if strand == "fw":
        ext = "/1"
    else:
        ext = "/2"
    fastq = ""
    for rownr in range(0, len(fastqDF.index)):
        fastq += str(fastqDF.index.values[rownr])+ext+"\n"
        fastq += str(fastqDF.iloc[rownr]["original_"+strand+"_seq"])+"\n"
        fastq += "+\n"
        fastq += str(fastqDF.iloc[rownr]["original_"+strand+"_seq_score"])+"\n"

    path = "data/"+sessionID+"/"+strand+"_file.fastq"
    with open(path,"w") as file:
        file.writelines(fastq)
    return path


def splitTSV(tsv_file):
    """Function to split the tsv file into options and data.

    Args:
        tsv_file: String of path to TSV file.

    Returns:
        options: String containing options.
        data: String containing data.
    
    """
    with open(tsv_file,"r") as file:
        lines = file.readlines()
        options = lines[0]
        data = "".join(lines[2:])

    return options, data


def configureOptions(options, fastqDF):
    """Update the FastQDataframe object with the options given.

    Args:
        options: A string with options extracted from tsv file.
        fastqDF: A FastQDataframe object that needs to be updated.

    """
    preprocess, blast = makeOptionsDict(options)
    if blast:
        fastqDF.set_blastOptions(blast)
    fastqDF.merge_preprocessOptionsDicts(preprocess)


def makeOptionsDict(options):
    """Convert options which were extracted from a tsv file to dictionaries.

    Args:
        options: A string containing options.

    Returns:
        preprocess_dict: A dictionary with options for preprocessing.
        blast_dict: A dictionary with options for BLAST.

    """
    preprocess_dict = {}
    blast_dict = {}
    for option in options.split("|"):
        if "=" in option:
            keyValue = option.split("=")
            try:
                value = int(keyValue[1].strip())
            except:
                raise ValueError
            key = keyValue[0].strip()
            if key.startswith("#"):
                key = key[1:]
            if key in ["minIdent", "reward", "penalty", "gapOpen", "gapExtend", "use_options"]:
                blast_dict[key] = value
            else:
                preprocess_dict[key] = value

    return preprocess_dict, blast_dict


def addMetaData(uuid, metaData):
    """Function to add the metadata to the fastqDataFrame object and check headers.

    Args:
        metaData: = pandas dataframe object
    
    Returns:
        True: If the metadata contains columns added by blast or sunburst.
        False: If the metadata does not contain columns added by blast or sunburst.

    """
    fastq_df = FastQDataframe.load_pickle("data/" + uuid + "/pickle.pkl")
    fastq_df.combine_dataframes(metaData)
    fastq_df.to_pickle(path='data/'+uuid+'/', filename='pickle')
    if "fw_full_tax" in metaData.columns or "rv_full_tax" in metaData.columns:
        return True
    elif "fw_visual_flag" in metaData.columns:
        raise IndexError('Error: TSV contains visual_flag column without Blast columns.\n'+
                         'Please check the TSV file or start a new session.')
    else:
        return False


def setFlags(flags, fastq_df):
    """Function to put the flags back from an old session.

    """
    df = fastq_df.get_dataframe()
    for column in flags.columns:
        df[column] = flags[column]
