import pandas as pd
import numpy as np
import json
import uuid
import pickle
from os.path import abspath
from app.core.objects.FastQData import FastQData
from app.core.preprocessing.calculations import get_sequence_length


class FastQDataframe:

    __df = pd.DataFrame
    __uuid = ""

    # Dictionaries to save the used options in:
    __preprocOptions = {}
    __originalPreprocOptions = {}
    __blastOptions = {}

    def __init__(self, fw_data: FastQData=None, rvc_data: FastQData=None, fw_items: list=None, fw_indices: list=None,
                 rvc_items: list=None, rvc_indices: list=None, df_id: str=None, fw_dict=None, rv_dict=None, from_blast=None):
        """Creates a FastQ dataframe.

        A FastQDataframe can be constructed multiple ways. The most straight forward way is to use FastQData objects,
        But it is also possible to use list items.
    
        Args:
            fw_data: FastQData.
            rvc_data: FastQData.
            fw_items: List with forward sequence and quality score items.
            fw_indices: List with identifiers of the forward items (fastQ headers).
            rvc_items: List with reverse complement sequence and quality score items.
            rvc_indices: List with identifiers of the reverse complement items (fastQ headers).
            df_id: Unique identifier for a FastQDataframe instance. Is used for writing this object to a pickle.
            fw_dict: ???
            rv_dict: ???
            from_blast: ???

        """
        if df_id:
            self.__uuid = df_id
        else:
            self.__uuid = str(uuid.uuid1())
        if fw_items and fw_indices and rvc_items and rvc_indices:
            self.__from_lists(fw_items, fw_indices, rvc_items, rvc_indices)
        elif fw_data and rvc_data:
            self.__from_fastqdata(fw_data, rvc_data)
        elif fw_dict and rv_dict:   
            self.__from_dict(fw_dict, rv_dict)
        elif from_blast:
            self.__from_blast()            
        else:
            raise Exception("No input data given, please use either FastQData objects or list objects")
            
    def __from_blast(self):
        fw_rvc_columns = []
        for prefix in ['fw_', 'rv_']:
            fw_rvc_columns.extend([prefix+'accession', prefix+'full_tax', prefix+'id', prefix+'bit', prefix+'coverage_length',
                       prefix+'Kingdom', prefix+'Phylum', prefix+'Class', prefix+'Order', prefix+'Family', prefix+'Genus', prefix+'Species', prefix + 'Strain'])
        fw_rvc_columns.extend(['paired_flag', 'fw_a_perc_flag', 'fw_t_perc_flag','fw_g_perc_flag', 'fw_c_perc_flag',
                   'rv_a_perc_flag', 'rv_t_perc_flag', 'rv_g_perc_flag', 'rv_c_perc_flag', 'fw_seq_len_flag',
                   'rv_seq_len_flag', 'identity_flag', 'fw_seq', 'rvc_seq' ])
        self.__df = pd.DataFrame(columns=fw_rvc_columns)

    # Cannot find what this function is used for.
    @DeprecationWarning
    def __from_lists(self, items_fw: list, indices_fw: list, items_rvc: list, indices_rvc: list):
        """Creates a pandas dataframe from list items.
        
        Args:
            items_fw: ?
            indices_fw: ?
            items_rvc: ?
            indices_rvc: ?

        """
        fw_columns = ['fw_seq', 'fw_seq_score']
        fw_df = pd.DataFrame(data=items_fw, columns=fw_columns, index=indices_fw)

        rvc_columns = ['rvc_seq', 'rvc_seq_score']
        rvc_df = pd.DataFrame(data=items_rvc, columns=rvc_columns, index=indices_rvc)

        self.__df = fw_df.join(rvc_df)

    def __from_fastqdata(self, fw_data: FastQData, rvc_data: FastQData):
        """Create a pandas dataframe from FastQData objects.
        
        Args:
            fw_data: ?
            rvc_data: ?

        """
        fw_columns = ['index_sequence','original_fw_seq', 'original_fw_seq_score',
                      'fw_seq',    'fw_seq_score', 'fw_seq_length', 
                      'fw_A_perc', 'fw_T_perc',    'fw_G_perc', 
                      'fw_C_perc', 'fw_Q-score']
        fw_df = pd.DataFrame(columns = fw_columns, 
                             index   = fw_data.get_index_list())

        if len(fw_data.get_info_index_list()) > 0:
            fw_df['index_sequence'] = fw_data.get_info_index_list()
            fw_df["index_sequence"] = fw_df["index_sequence"].str.split(":").str[3]
        else:
            fw_df.drop(['index_sequence'], axis=1, inplace=True)
        fw_df['fw_seq'] = fw_data.get_sequence_list()
        fw_df['fw_seq_score'] = fw_data.get_q_score_list()
        fw_df['original_fw_seq'] = fw_data.get_sequence_list()
        fw_df['original_fw_seq_score'] = fw_data.get_q_score_list()
        
        rvc_columns = ['original_rv_seq', 'original_rvc_seq', 'original_rv_seq_score',
                      'rv_seq',       'rvc_seq',   'rv_seq_score', 
                      'rv_seq_length', 'rv_A_perc', 'rv_T_perc', 
                      'rv_G_perc',     'rv_C_perc', 'rv_Q-score']
        rvc_df = pd.DataFrame(columns = rvc_columns, 
                              index   = rvc_data.get_index_list())
        
        rvc_df['rv_seq'] = rvc_data.get_sequence_list()
        rvc_df['rvc_seq'] = rvc_data.get_complement_list()
        rvc_df['rv_seq_score'] = rvc_data.get_q_score_list()
        rvc_df['original_rv_seq'] = rvc_data.get_sequence_list()
        rvc_df['original_rvc_seq'] = rvc_data.get_complement_list()
        rvc_df['original_rv_seq_score'] = rvc_data.get_q_score_list()

        joined_df = fw_df.join(rvc_df)
        
        # TODO: Check if new columns need to be added to the table for future features.
        for newColumn in ['flagged','fw_flagged','rv_flagged','paired_flag',
                        'fw_a_perc_flag','fw_t_perc_flag', 'fw_g_perc_flag',  
                        'fw_c_perc_flag','rv_a_perc_flag', 'rv_t_perc_flag',  
                        'rv_g_perc_flag','rv_c_perc_flag', 'fw_seq_len_flag', 
                        'rv_seq_len_flag','identity_flag', 'fw_quality_flag', 
                        'rv_quality_flag','original_contains_unknown_symbols']:
            # Fill the rows of the new column with False.
            joined_df[newColumn] = False  
        
        self.__df = joined_df
        
        self.__df['original_contains_unknown_symbols'] = list(map(self.__has_unknown, self.__df[["original_fw_seq","original_rv_seq"]].values.tolist()))
        
        self.calcQuality("fw")
        self.calcQuality("rv")

    @staticmethod    
    def __has_unknown(row):
        """Checks if the original seq had unknown symbols (these were replaced with N in previous processes)."""
        def check_N(item):
            if not pd.isnull(item):
                if "N" in item.upper():
                    return True
            return False
        
        for item in row:
            if check_N(item):
                return True
        return False

    def __from_dict(self, fw_data, rv_data):
        """Create a pandas dataframe from FastQData objects.

        Args:
            fw_data: ?
            rvc_data: ?

        """
        fw_columns = ['fw_seq',    'fw_seq_score', 'fw_seq_length', 
                      'fw_A_perc', 'fw_T_perc',    'fw_G_perc', 
                      'fw_C_perc', 'fw_Q-score']
        fw_df = pd.DataFrame(columns=fw_columns, index=fw_data['index'])

        fw_df['fw_seq'] = fw_data['seq']
        fw_df['fw_seq_score'] = fw_data['qual']


        rv_columns = ['rv_seq',        'rvc_seq',   'rv_seq_score', 
                      'rv_seq_length', 'rv_A_perc', 'rv_T_perc', 
                      'rv_G_perc',     'rv_C_perc', 'rv_Q-score']
        rv_df = pd.DataFrame(columns=rv_columns, index=rv_data['index'])
        rv_df['rv_seq'] = rv_data['seq']
        rv_df['rvc_seq'] = rv_data['comp']
        rv_df['rv_seq_score'] = rv_data['qual']

        del fw_data
        del rv_data

        joined_df = fw_df.join(rv_df)

        del fw_df
        del rv_df

        for newColumn in ['flagged',        'paired_flag',     'fw_a_perc_flag',
                          'fw_t_perc_flag', 'fw_g_perc_flag',  'fw_c_perc_flag', 'fw_Q-score',
                          'rv_a_perc_flag', 'rv_t_perc_flag',  'rv_g_perc_flag',
                          'rv_c_perc_flag', 'rv_Q-score', 'fw_seq_len_flag', 'rv_seq_len_flag',
                          'identity_flag']:
            joined_df[newColumn] = False  # Fill the rows of the new column with False.
       
        self.__df = joined_df
        self.calcQuality("fw")
        self.calcQuality("rv")

    def set_df(self, df: pd.DataFrame):
        self.__df = df

    def update_by_list(self, data: list, columns: list):
        # Items of data must be the same length as the number of indexes in the dataframe.
        for i, column in enumerate(columns):
            self.__df[column] = [data][i]

    def update_by_series(self, data: pd.Series):
        self.__df.update(data)

    def update_by_dataframe(self, data: pd.DataFrame):
        self.__df.update(data)

    def get_dataframe(self):
        return self.__df

    def combine_dataframes(self, df: pd.DataFrame):
        self.__df = pd.concat([self.__df, df], axis = 1)

    def get_list_columns(self):
        return self.__df.columns

    def get_column(self, column_name: str):
        return self.__df[column_name]

    def get_columns(self, column_names: list):
        return self.__df[column_names]

    def get_uuid(self):
        return self.__uuid

    def filter_between(self, value1, value2, columns: list):
        filtered_df = self.__df.copy()
        for column in columns:
            filtered_df = filtered_df[filtered_df[column].between(value1, value2, inclusive=True)]
        return filtered_df

    def filter_greater_than(self, value, columns: list):
        filtered_df = self.__df.copy()
        for column in columns:
            filtered_df = filtered_df[filtered_df[column] > value]
        return filtered_df

    def filter_smaller_than(self, value, columns: list):
        filtered_df = self.__df.copy()
        for column in columns:
            filtered_df = filtered_df[filtered_df[column] < value]
        return filtered_df

    def filter_equals(self, value, columns: list):
        filtered_df = self.__df.copy()
        for column in columns:
            filtered_df = filtered_df[filtered_df[column] == value]
        return filtered_df

    @staticmethod
    def is_filtered(row):
        try:
            if row['flagged'] == False:
                return True
        except KeyError:
            return False  

    def flag_any(self, fwColumns: list=['paired_flag','fw_seq_len_flag','fw_a_perc_flag', 'fw_t_perc_flag','fw_g_perc_flag', 'fw_c_perc_flag','fw_quality_flag', 'identity_flag'],
                        rvColumns: list=['paired_flag','rv_seq_len_flag','rv_a_perc_flag', 'rv_t_perc_flag','rv_g_perc_flag', 'rv_c_perc_flag','rv_quality_flag', 'identity_flag'], 
                        flaggedColumns : list=['fw_flagged', 'rv_flagged'], axis=1,
                        fw_flag_column = 'fw_flagged', rv_flag_column = 'rv_flagged', flagged = 'flagged'):
        self.__df[fw_flag_column] = self.__df[fwColumns].any(axis)
        self.__df[rv_flag_column] = self.__df[rvColumns].any(axis)
        self.__df[flagged] = self.__df[flaggedColumns].any(axis)

    def flag_between(self, filter_on, value1, value2, column='flagged'):
        self.__df[column] = self.__df[filter_on].apply(lambda row: False if value1 <= row <= value2 else True)

    @staticmethod
    def __do_between(row, column, value1, value2):
        if value1 <= row[column] <= value2:
            return False
        else:
            return True

    def flag_smaller_than(self, filter_on, value, column='flagged'):
        self.__df[column] = self.__df[filter_on].apply(lambda row: False if row <= value else True)

    @staticmethod
    def __do_smaller(row, column, value):
        if row[column] <= value:
            return False
        else:
            return True

    def flag_greater_than(self, filter_on, value, column='flagged'):
        self.__df[column] = self.__df[filter_on].apply(lambda row: False if row >= value else True)

    @staticmethod
    def __do_greater(row, column, value):
        if row[column] >= value:
            return False
        else:
            return True

    def flag_reset_column(self, column, value = False):
        self.__df[column] = value

    def flag_equals(self, filter_on, value, column='flagged'):
        self.__df[column] = self.__df[filter_on].apply(lambda row: False if row == value else True)

    def flag_indexes(self,indexList, column = 'flagged'):
        self.__df.loc[indexList,column] = True 

    @staticmethod
    def __do_equals(row, column, value):
        if row[column] == value:
            return False
        else:
            return True

    def to_json(self):
        return self.__df.to_json()

    def to_csv(self, path):
        return self.__df.to_csv(path)

    def to_pickle(self, path='data/', filename=None):
        if filename:
            with open(abspath(path+filename+".pkl"), "wb") as f:
                pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        else:
            with open(abspath(path+self.__uuid+".pkl"), "wb") as f:
                pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load_pickle(file):
        with open(abspath(file), 'rb') as f:
            df = pickle.load(f)
            if type(df) is FastQDataframe:
                return df
            else:
                raise TypeError('Pickle is not a FastQDataframe')

    def get_preprocessOptions(self):
        return self.__preprocOptions

    def get_originalPreprocessOptionsJSON(self):
        return self.__originalPreprocOptions

    def get_blastOptions(self):
        return self.__blastOptions

    #TODO see what options to delete later on.
    def set_preprocessOptions(self, paired_read_percentages = 0,
                                 min_A_perc  = 0,    min_C_perc     = 0,
                                 min_G_perc     = 0,   min_T_perc  = 0,    max_A_perc     = 100,
                                 max_C_perc     = 100, max_G_perc  = 100,  max_T_perc     = 100):  

        min_A_perc = int(np.floor(self.__df[['fw_A_perc','rv_A_perc']].min(axis=1).min()))
        min_T_perc = int(np.floor(self.__df[['fw_T_perc','rv_T_perc']].min(axis=1).min()))
        min_G_perc = int(np.floor(self.__df[['fw_G_perc','rv_G_perc']].min(axis=1).min()))
        min_C_perc = int(np.floor(self.__df[['fw_C_perc','rv_C_perc']].min(axis=1).min()))
        max_A_perc = int(np.ceil(self.__df[['fw_A_perc','rv_A_perc']].max(axis=1).max()))
        max_T_perc = int(np.ceil(self.__df[['fw_T_perc','rv_T_perc']].max(axis=1).max()))
        max_G_perc = int(np.ceil(self.__df[['fw_G_perc','rv_G_perc']].max(axis=1).max()))
        max_C_perc = int(np.ceil(self.__df[['fw_C_perc','rv_C_perc']].max(axis=1).max()))

        self.__preprocOptions = {"paired_read_percentages" : paired_read_percentages,
                                 "min_A_perc"  : min_A_perc,  "min_C_perc" : min_C_perc,
                                 "min_G_perc"   : min_G_perc,   "min_T_perc"  : min_T_perc,  "max_A_perc" : max_A_perc,
                                 "max_C_perc"   : max_C_perc,   "max_G_perc"  : max_G_perc,  "max_T_perc" : max_T_perc}
        self.__originalPreprocOptions = self.__preprocOptions.copy()
    
    def getForwardStats(self, percentageDict):
        self.__df = get_sequence_length(self.__df, True)

        qualityCutoff = self.__df['fw_Q-score'] >= percentageDict['fw_Q_cutoff']

        minA = self.__df['fw_A_perc'] >= percentageDict['min_A_perc']
        minT = self.__df['fw_T_perc'] >= percentageDict['min_T_perc']
        minG = self.__df['fw_G_perc'] >= percentageDict['min_G_perc']
        minC = self.__df['fw_C_perc'] >= percentageDict['min_C_perc']

        maxA = self.__df['fw_A_perc'] <= percentageDict['max_A_perc']
        maxT = self.__df['fw_T_perc'] <= percentageDict['max_T_perc']
        maxG = self.__df['fw_G_perc'] <= percentageDict['max_G_perc']
        maxC = self.__df['fw_C_perc'] <= percentageDict['max_C_perc']
        df = self.__df[qualityCutoff&minA&minT&minG&minC&maxA&maxT&maxG&maxC]
        fwTotal = df['fw_flagged'][df.fw_flagged == False].count()
        avgLength = self.__df.loc[self.__df['fw_flagged'] == False, 'fw_seq_length'].mean()

        return fwTotal, avgLength
    
    def getReverseStats(self, percentageDict):
        qualityCutoff = self.__df['rv_Q-score'] >= percentageDict['rv_Q_cutoff']

        minA = self.__df['rv_A_perc'] >= percentageDict['min_A_perc']
        minT = self.__df['rv_T_perc'] >= percentageDict['min_T_perc']
        minG = self.__df['rv_G_perc'] >= percentageDict['min_G_perc']
        minC = self.__df['rv_C_perc'] >= percentageDict['min_C_perc']

        maxA = self.__df['rv_A_perc'] <= percentageDict['max_A_perc']
        maxT = self.__df['rv_T_perc'] <= percentageDict['max_T_perc']
        maxG = self.__df['rv_G_perc'] <= percentageDict['max_G_perc']
        maxC = self.__df['rv_C_perc'] <= percentageDict['max_C_perc']
        df = self.__df[qualityCutoff&minA&minT&minG&minC&maxA&maxT&maxG&maxC]
        rvTotal = df['rv_flagged'][df.rv_flagged == False].count()
        avgLength = self.__df.loc[self.__df['rv_flagged'] == False, 'rv_seq_length'].mean()
        return rvTotal, avgLength

    def getMaxLength(self, strand = ""):
        maxFw = max(self.__df[self.__df["fw_seq_length"].notnull()]["fw_seq_length"])
        # maxFw = self.__df.original_fw_seq.map(len).max()
        if strand == "fw":
            return maxFw
        maxRv = max(self.__df[self.__df["rv_seq_length"].notnull()]["rv_seq_length"])
        # maxRv = self.__df.original_rv_seq.map(len).max()
        if strand == "rv":
            return maxRv
        return max(maxFw,maxRv)

    def calcQuality(self,strand):
        self.__df[strand+"_Q-score"] = list(map(lambda qual: 
                                                sum(list(map(lambda s: 
                                                                ord(s) - 33, 
                                                            [i for i in qual])
                                                        )
                                                    ) / len(qual)  if not pd.isnull(qual) and len(qual) != 0 else 0, 
                                              list(self.__df.loc[:,strand+"_seq_score"])
                                              )
                                        )
    

    def get_read_length(self, strand = "fw"):
        return list(map(lambda x: x if not pd.isnull(x) else 0, self.__df[strand+"_seq_length"]))

    def getMinLength(self, strand = ""):
        minFw = min(self.__df[self.__df["fw_seq_length"].notnull()]["fw_seq_length"])
        if strand == "fw":
            return minFw
        minRv = min(self.__df[self.__df["rv_seq_length"].notnull()]["rv_seq_length"])
        if strand == "rv":
            return minRv
        return min(minFw,minRv)

    def getAvgQuality(self):
        avgFw = round(self.__df.loc[self.__df['fw_flagged'] == False, 'fw_Q-score'].mean(), 2)
        avgRv = round(self.__df.loc[self.__df['rv_flagged'] == False, 'rv_Q-score'].mean(), 2)
        return avgFw, avgRv

    def getMaxQuality(self):
        fwMaxQuality = np.floor(self.__df['fw_Q-score'].max())
        rvMaxQuality = np.floor(self.__df['rv_Q-score'].max())
        return fwMaxQuality, rvMaxQuality

    def flag_qual(self, cutoff, strand):
        self.__df[strand+'_quality_flag'] = (self.__df[strand+'_Q-score'] <= cutoff)

    def flag_test(self, preDict):
        self.__df['fw_quality_flag'] = (self.__df['fw_Q-score'] <= preDict['fw_Q_cutoff'])
        self.__df['rv_quality_flag'] = (self.__df['rv_Q-score'] <= preDict['rv_Q_cutoff'])

    def getTotalReads(self):
        fwTotal = self.__df['fw_flagged'][self.__df.fw_flagged == False].count()
        rvTotal = self.__df['rv_flagged'][self.__df.rv_flagged == False].count()
        return (fwTotal, rvTotal)
    
    def getUnflaggedReads(self):
        return self.__df['flagged'][self.__df.flagged == False].count()
    
    def setRangeOptions(self, rangeDict):
        self.__preprocOptions.update(rangeDict)
    
    def set_blastOptions(self, blastOptions):
        # TODO: If you add, change or remove options please also do this in app/containr/preprocessing/parser.py @ makeOptionsDict()
        self.__blastOptions = {"minIdent" : blastOptions["minIdent"], "reward": blastOptions["reward"], 
                            "penalty": blastOptions["penalty"],  "gapOpen": blastOptions["gapOpen"],  "gapExtend": blastOptions["gapExtend"]}


    def merge_blastOptionsDicts(self,blastDict):
        self.__blastOptions.update(blastDict)
    
    def merge_preprocessOptionsDicts(self,preDict):
        self.__preprocOptions.update(preDict)
