import io
import json
import logging
import os

from collections import Counter

from flask import jsonify, request, send_file, abort

import pandas as pd

from app.api import bp

from app.core.containr.Json_parser import receive_json_ids
from app.core.containr.index import containr_main
from app.core.containr.parse_blast_results import *
from app.core.containr.visualise_identity import *
from app.core.objects.FastQDataframe import FastQDataframe
from app.core.preprocessing.calculations import make_alignment, get_sequence_length
from app.core.utils.to_json import seq_length_to_json, perc_count_to_json, get_paired_percentage_to_json, \
    nucleotide_percentages_to_json


@bp.route('/sequenceLength', methods=['POST'])
def sequence_length():
    """
    Endpoint for getting filtered sequence length data from the FastQDataframe
    :return: json object containing forward and reverse sequence length data
    """
    try:
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        min_seq_len = int(request.form.get('min_seq_len'))
        max_seq_len = int(request.form.get('max_seq_len'))
        fastq_df.flag_between("fw_seq_length", min_seq_len, max_seq_len, 'fw_seq_len_flag')
        fastq_df.flag_between("rv_seq_length", min_seq_len, max_seq_len, 'rv_seq_len_flag')
        fastq_df.flag_any()
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')
        subset = fastq_df.filter_equals(False, ['flagged'])
        response = seq_length_to_json(subset[['fw_seq_length', 'rv_seq_length']])
        return response
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/nucleotide', methods=['POST'])
def nucleotide():
    """
    Endpoint for getting filtered data about nucleotide ratios
    Added: Now possible to filter on a specific nucleotide
    :return: json object
    """
    try:
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        min_A_perc = int(request.form.get('minAValue'))
        min_T_perc = int(request.form.get('minTValue'))
        min_G_perc = int(request.form.get('minGValue'))
        min_C_perc = int(request.form.get('minCValue'))
        max_A_perc = int(request.form.get('maxAValue'))
        max_T_perc = int(request.form.get('maxTValue'))
        max_G_perc = int(request.form.get('maxGValue'))
        max_C_perc = int(request.form.get('maxCValue'))

        fastq_df.flag_between('fw_A_perc', min_A_perc, max_A_perc, 'fw_a_perc_flag')
        fastq_df.flag_between('fw_T_perc', min_T_perc, max_T_perc, 'fw_t_perc_flag')
        fastq_df.flag_between('fw_G_perc', min_G_perc, max_G_perc, 'fw_g_perc_flag')
        fastq_df.flag_between('fw_C_perc', min_C_perc, max_C_perc, 'fw_c_perc_flag')

        fastq_df.flag_between('rv_A_perc', min_A_perc, max_A_perc, 'rv_a_perc_flag')
        fastq_df.flag_between('rv_T_perc', min_T_perc, max_T_perc, 'rv_t_perc_flag')
        fastq_df.flag_between('rv_G_perc', min_G_perc, max_G_perc, 'rv_c_perc_flag')
        fastq_df.flag_between('rv_C_perc', min_C_perc, max_C_perc, 'rv_g_perc_flag')
        fastq_df.flag_any()
        subset = fastq_df.filter_equals(False, ['flagged'])
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')

        fw_json = nucleotide_percentages_to_json(subset[['fw_A_perc', 'fw_T_perc', 'fw_C_perc', 'fw_G_perc']])
        rv_json = nucleotide_percentages_to_json(subset[['rv_A_perc', 'rv_T_perc', 'rv_C_perc', 'rv_G_perc']])

        response = json.dumps({"fw_json": fw_json, "rvc_json": rv_json})
        return response
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except IndexError as e:
        logging.exception(e)
        abort(411)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/nucleotidePercentage', methods=['POST'])
def nucleotidePercentage():
    try:
        session_id = session['id']      # used to know in which directory the data can be found
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl") # main object with the dataframe
        df = fastq_df.get_dataframe()  # dataframe containing the data

        # Used to determine which strand should be used
        if request.form.get("forwardRead") == "true":
            strand = "fw"
        else:
            strand = "rv"

        #Cut off reads with quality cutoff
        quality_cutoff = int(request.form.get("qualityCutoff"))
        fastq_df.flag_qual(quality_cutoff,strand)
        fastq_df.flag_any()

        df = df[df["flagged"] == False]

        dictRange = {}

        # Ranges obtained from the GUI
        min_read_range = int(request.form.get("minReadRange"))
        max_read_range = int(request.form.get("maxReadRange"))

        dictRange["minRange" + strand] = min_read_range
        dictRange["maxRange" + strand] = max_read_range
        fastq_df.setRangeOptions(dictRange)

        return_value = {"A": [], "C": [], "G": [], "T": [],"max": 0,"average":[],"range":[]} # dict with counts of each nucleotide and quality per position

        df = df[df["flagged"] == False]
        
        # Ranges obtained from the GUI
        min_read_range = int(request.form.get("minReadRange"))
        max_read_range = int(request.form.get("maxReadRange"))

        return_value = {"A": [], "C": [], "G": [], "T": [],"max": 0,"average":[],"range":[]} # dict with counts of each nucleotide and quality per position
        
        # Make a dictionary for each position with the count of te nucleotide
        for index in range(min_read_range - 1, max_read_range):
            try:
                countDict = Counter(df.loc[:,"original_"+strand+"_seq"].str[index])
                # Check if the max value needs to be updated
                maxCount = max(countDict.values())
                if maxCount > return_value["max"]:
                    return_value["max"] = maxCount
                
                # Fill dict if it miss values
                for key in set(["A","C","G","T"]).difference(set(countDict.keys())):
                    countDict[key] = 0
            except IndexError:
                countDict = {"A":0,"C":0,"G":0,"T":0}
            
            # Write values to value dict
            return_value["A"].append({"x": index + 1, "y": countDict["A"]})  # x = position, y = counter
            return_value["C"].append({"x": index + 1, "y": countDict["C"]})
            return_value["G"].append({"x": index + 1, "y": countDict["G"]})
            return_value["T"].append({"x": index + 1, "y": countDict["T"]})

            # Determine the quality scores average and mean
            try:
                qualList = list(map(lambda qual: ord(qual) - 33 if not pd.isnull(qual) else 0, df.loc[:,"original_"+strand+"_seq_score"].str[index]))
            except IndexError:
                qualList = [0]

            return_value["average"].append({"x": index + 1, "y": sum(qualList) / len(qualList)})
            return_value["range"].append({"x": index + 1, "y": [min(qualList), max(qualList)]})

    	
        # Adjust the dataframe so the other graphs can use it
        df = fastq_df.get_dataframe()
        df[strand+"_seq"] = df.loc[:,"original_"+strand+"_seq"].str[min_read_range - 1 : max_read_range]
        df[strand+"_seq_score"] = df.loc[:,"original_"+strand+"_seq_score"].str[min_read_range - 1 : max_read_range]
        if strand == "rv":
            df[strand+"c_seq"] = df.loc[:,"original_"+strand+"c_seq"].str[-max_read_range : - (min_read_range - 1)]
        # df = get_sequence_length(df,True)
        fastq_df.set_df(df)
        fastq_df.calcQuality("fw")
        fastq_df.calcQuality("rv")
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')

        return jsonify(return_value)
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)

@bp.route('/update_stats', methods=['POST'])
def update_stats():
    """
    Endpoint for update the statistics on the preprocess page
    :return: json object
    """
    
    try:
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        optDict = { 'fw_Q_cutoff' : int(request.form.get('fwQCut')), 'rv_Q_cutoff' : int(request.form.get('rvQCut')),
                    'min_A_perc' : int(request.form.get('minAValue')), 'min_T_perc' : int(request.form.get('minTValue')),
                    'min_G_perc' : int(request.form.get('minGValue')), 'min_C_perc' : int(request.form.get('minCValue')),
                    'max_A_perc' : int(request.form.get('maxAValue')), 'max_T_perc' : int(request.form.get('maxTValue')),
                    'max_G_perc' : int(request.form.get('maxGValue')), 'max_C_perc' : int(request.form.get('maxCValue')),
                    'minSeq' : int(request.form.get('minSeqLen')), 'maxSeq' : int(request.form.get('maxSeqLen'))}
      
        fastq_df.flag_test(optDict)
        fastq_df.flag_any()
        
        fwCount, fwAvg = fastq_df.getForwardStats(optDict)
        rvCount, rvAvg = fastq_df.getReverseStats(optDict)
        avgQFw, avgQRv = fastq_df.getAvgQuality()
        unflaggedReads = fastq_df.getUnflaggedReads()

        response = json.dumps({"fwReads": str(fwCount), "rvReads": str(rvCount), 
                                "avgQualityFw": str(avgQFw), "avgQualityRv": str(avgQRv), "avgLengthFw": str(fwAvg), 
                                "avgLengthRv": str(rvAvg), "unflaggedReads": str(unflaggedReads)})
         
        fastq_df.merge_preprocessOptionsDicts(optDict)
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')

        return response
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)

@bp.route('/quality', methods=['POST'])
def quality():
    """
    Endpoint for getting a subset of the dataframe filtered on a quality score
    :return: json object
    """
    try:
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        response = request.form.get('qValue')
        return jsonify(response)
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/paired', methods=['POST'])
def paired():
    """
    Endpoint for getting a subset containing only paired reads
    :return: json object for a image
    """
    try:
        session_id = session['id']
        filter_paired = True if request.form.get('FilterPaired') == "true" else False
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        if filter_paired:
            fastq_df.flag_equals('paired', filter_paired, 'paired_flag')
        else:
            fastq_df.flag_reset_column("paired_flag")
        fastq_df.flag_any()
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')
        response = get_paired_percentage_to_json(fastq_df.get_dataframe())
        return response
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/identity' ,methods=['POST'])
def identity():
    """
    Endpoint for filtering the dataframe on identity score and getting json data for an image
    :return: json object
    """
    
    try:
        # get parameters from gui
        iden_perc = int(request.form.get('paired_read_minimum_percentage'))
        matchScore = int(request.form.get("paired_read_match"))
        mismatchScore = int(request.form.get("paired_read_mismatch"))
        open_gap_score = int(request.form.get("paired_read_opening_gap"))
        extend_gap_score = int(request.form.get("paired_read_extend_gap"))

        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")

        fastq_df.set_df(make_alignment(fastq_df.get_dataframe(), matchScore, mismatchScore, open_gap_score, extend_gap_score))

        fastq_df.flag_greater_than('overlap_identity_perc', iden_perc, 'identity_flag')
        fastq_df.flag_any()
        filtered_df = fastq_df.filter_equals(False, ['flagged'])
        subset = filtered_df[['overlap_identity_perc']].astype(float).round().groupby(['overlap_identity_perc']).overlap_identity_perc.count()
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')
        response = perc_count_to_json(subset)
        
        return jsonify(response)
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/export_tsv', methods=['GET'])
def export_tsv():
    """
    Endpoint for returning a TSV file of the dataframe
    :return: tsv file
    """

    try:
        session_id = session['id']
        
        fastq_df = FastQDataframe.load_pickle("data/"+ session_id + "/pickle.pkl")
        df = fastq_df.get_dataframe()

        buffer = io.StringIO()
        optionsText = "#"
        for options in [fastq_df.get_preprocessOptions(),fastq_df.get_blastOptions()]:
            for key in options.keys():
                optionsText += key+"="+str(options[key])+"|"

        buffer.write(optionsText
                     + "\n#columns with flag in the name; True means it's filtered out, False means it's included in future analyses. Column Original_contains_unknown_symbols stands for the fact that if the original sequences contained symbols that did not match A, C, G or T. These unknown symbols are replaced with the universal symbol N. The original sequence columns already underwent this change to safe time for the application.\n")
        
        df.to_csv(buffer, sep='\t', index=True, header=True)
        
        mem = io.BytesIO()
        mem.write(buffer.getvalue().encode('utf-8'))
        mem.seek(0)
        buffer.close()

        return send_file(mem,
                         mimetype='text/tab-separated-values',
                         attachment_filename='ConTaInR_export_data.tsv',
                         as_attachment=True,
                         cache_timeout=-1)
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/blastDuration', methods=['POST'])
def blastDuration():
    try:
        use_options = request.form.get("use_options") == "true"
        session_id = session['id']
        fw_size = os.path.getsize(f"data/{session_id}/fw_file.fastq")
        rv_size = os.path.getsize(f"data/{session_id}/rv_file.fastq")
        # Some complicated calculations led to the following multiplication value but feel free to edit as the speed of the application changes.
        duration = 0.0000044 * (fw_size + rv_size)
        if use_options:
            duration *= 2.8
        response = json.dumps({"duration": round(duration)})

        return response
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)


@bp.route('/BLAST', methods=['POST'])
def BLAST():
    """
    Endpoint for executing BLAST
    :return:
    """
    try:
        session_id = session['id']
        use_options = request.form.get("use_options") == "true"
        minIdent_value = int(request.form.get("minIdent"))
        reward_value = int(request.form.get("reward"))
        penalty_value = int(request.form.get("penalty")) * -1
        open_value = int(request.form.get("open"))
        extend_value = int(request.form.get("extend"))

        blastOptions = {"minIdent": minIdent_value, "reward": reward_value, "penalty": penalty_value, "gapOpen": open_value, "gapExtend": extend_value, "use_options": use_options}

        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl") # Get object with all the data
        containr_main("app/core/containr/Ref_files/ref_file.txt", fastq_df, blastOptions)  

        fastq_df.set_blastOptions(blastOptions)
        main_df = fastq_df.get_dataframe()
        if "fw_visual_flag" not in main_df.columns:
            main_df["fw_visual_flag"] = False
            main_df["rv_visual_flag"] = False
                                     
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')
        return json.dumps({"success": True})
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)

@bp.route('/load_records', methods=['POST'])        
def load_records():
    """
    Endpoint for extending the dataframe with blast output
    :return:
    """
    try:
        session_id = session['id']
        update_dataframe = str(request.form.get('update_dataframe'))
        fastq_df = (FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl"))
        return json.dumps({"success": True})
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)
    

@bp.route('/sunburst', methods=['POST'])
def sunburst():
    """
    Endpoint for getting the data for the sunburst
    :return:
    """
    try:
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        df = fastq_df.get_dataframe()
        remake = True if request.form.get('remake') == "True" else False
        if remake:
            df["fw_visual_flag"] = False
            df["rv_visual_flag"] = False
            fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')
        response = {} 
        response['sun_json'] = convert_results(df)
        return json.dumps(response)
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)

@bp.route('/delete_node', methods=['POST'])
def delete_node():
    """
    Endpoint for deleting a segment of the sunburst
    :return:
    """
    try:
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        jsondata = str(request.form.get('data_json'))
        fw_index, rv_index = receive_json_ids(fastq_df.get_dataframe(),jsondata, just_headers=True)
        fastq_df.flag_indexes(fw_index,'fw_visual_flag') 
        fastq_df.flag_indexes(rv_index,'rv_visual_flag')
        fastq_df.to_pickle(path='data/' + session_id + '/', filename='pickle')
        return json.dumps({"success": True})
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)

@bp.route('/bacterial', methods=['POST'])
def bacterial_visualization():
    """
    Endpoint for getting data when a section of the sunburst is clicked
    :return: JSON object containing 3 JSON object for the 3 summary graphs
    """
    try:
        jsondata = request.form.get("data_json")
        session_id = session['id']
        fastq_df = FastQDataframe.load_pickle("data/" + session_id + "/pickle.pkl")
        response = receive_json_ids(fastq_df.get_dataframe(), jsondata)
        return response
    except KeyError as e:
        logging.exception(e)
        abort(400)
    except Exception as e:
        logging.exception(e)
        abort(500)
