
class FastQData:

    __index_list = []
    __index_info_list = []
    __sequence_list = []
    __complement_list = []
    __q_score_list = []
    __is_rev = False

    def __init__(self, index_list: list, info_index_list: list, sequence_list: list, complement_list: list , q_score_list: list, is_rev: bool=False):
        self.__index_list = index_list
        self.__index_info_list = info_index_list
        self.__sequence_list = sequence_list
        self.__complement_list = complement_list
        self.__q_score_list = q_score_list
        self.__is_rev = is_rev

    def get_index_list(self):
        return self.__index_list

    def get_info_index_list(self):
        return self.__index_info_list

    def get_sequence_list(self):
        return self.__sequence_list

    def get_complement_list(self):
        return self.__complement_list

    def get_q_score_list(self):
        return self.__q_score_list

    def get_fastq_data(self):
        return self.__to_dict()

    def __to_dict(self):
        fastq_data = {
            "Indices": self.__index_list,
            "Indices-info": self.__index_info_list,
            "Sequences": self.__sequence_list,
            "Q-scores": self.__q_score_list,
            "Is-reverse": self.__is_rev
        }
        return fastq_data

    def is_reverse(self):
        return self.__is_rev
