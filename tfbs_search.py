# coding=utf-8
from Bio import SeqIO
from math import log, exp
import math
import os
import subprocess

global_compila_main = ""
global_gre_pwm_source = ""
global_pseudocontador_seq = 0.01
global_gre_error = set()

def parse_freq_gen(genfreqfile):

    genfreq = open(genfreqfile, 'r')

    dict_freq_gen = {}
    for freq in genfreq:
        if freq[1:2] == "#":
            continue

        pos_virg1 = freq.find(",")
        pos_virg2 = freq.find(",", pos_virg1 + 1)
        pos_virg3 = freq.find(",", pos_virg2 + 1)
        pos_virg4 = freq.find(",", pos_virg3 + 1)

        freq_g = float(freq[pos_virg1 + 1:pos_virg2])
        freq_c = float(freq[pos_virg2 + 1:pos_virg3])
        freq_a = float(freq[pos_virg3 + 1:pos_virg4])
        freq_t = float(freq[pos_virg4 + 1:-1])

        dict_freq_bas = {}
        dict_freq_bas['A'] = freq_a
        dict_freq_bas['C'] = freq_c
        dict_freq_bas['G'] = freq_g
        dict_freq_bas['T'] = freq_t

        dict_freq_gen[freq[1:pos_virg1 - 1]] = dict_freq_bas

    genfreq.close()
    return dict_freq_gen


def parse_tf_par(list_tf_par, tf_parameters):
    dict_tf_par = {}
    list_tf_par_f = ["FT" + str(ft).zfill(7) for ft in list_tf_par]

    with open(tf_parameters, 'r') as h_tf_par:
        for tf_par_line in h_tf_par:
            data = tf_par_line.split("\t")
            if (data[1].rstrip() == ""):
                continue
            if (data[1] not in list_tf_par_f):
                continue
            dict_tf_par[data[1]] = {}
            dict_tf_par[data[1]]["TF_CODE"] = data[2].rstrip()
            dict_tf_par[data[1]]["TF_NAME"] = data[3].rstrip()
            dict_tf_par[data[1]]["TF_TECH"] = data[4].rstrip()
            dict_tf_par[data[1]]["TF_ORIG_NAME"] = data[5].rstrip()
            dict_tf_par[data[1]]["TF_TFBS_LEN"] = int(data[6].rstrip())
            dict_tf_par[data[1]]["TF_TFBS_LEN"] = int(data[6].rstrip())
            dict_tf_par[data[1]]["TF_SAMPLE_SIZE"] = int(data[7].rstrip())

    return dict_tf_par


def parse_pwm(list_exe_par, list_models):
    os.chdir(global_gre_pwm_source)
    list_models = [x for x in list_models if x != "N3"]
    dict_pwm= {}
    for tf_name in list_exe_par:
        dict_pwm[tf_name] = {}
        for model in list_models:

            dict_pwm[tf_name][model] = {}
            dict_pwm[tf_name][model]["PWM"] = {}

            pwm_filename = "FT" + str(tf_name).zfill(7) + "_" + model + ".pwm"
            with open(pwm_filename, 'r') as pwm_file:
                position = 0
                for pwm_line in pwm_file:
                    data = pwm_line.split("\t")
                    if (data[0] == "#"):
                        if (data[1] != "A"):
                            dict_pwm[tf_name][model][data[1][:-1]] = float(data[2].rstrip())
                        else:
                            dict_pwm[tf_name][model]["PWM"]["A"] = {}
                            dict_pwm[tf_name][model]["PWM"]["C"] = {}
                            dict_pwm[tf_name][model]["PWM"]["G"] = {}
                            dict_pwm[tf_name][model]["PWM"]["T"] = {}
                    elif (data[0] == "*"):
                        dict_pwm[tf_name][model]["PWM"]["A"][position] = float(data[1].rstrip())
                        dict_pwm[tf_name][model]["PWM"]["C"][position] = float(data[2].rstrip())
                        dict_pwm[tf_name][model]["PWM"]["G"][position] = float(data[3].rstrip())
                        dict_pwm[tf_name][model]["PWM"]["T"][position] = float(data[4].rstrip())
                        position += 1
    return dict_pwm


def parse_gre(list_exe_par, gre_parameters):
    dict_gre= {}
    for ft in list_exe_par:
        tf_name = "FT" + str(ft).zfill(7)
        ft_found = False
        with open(gre_parameters, 'r') as gre_par:
            for gre_line in gre_par:
                if (gre_line.find(tf_name) < 0) and not ft_found:
                    continue
                ft_found = True
                data = gre_line.split("\t")
                if (data[0] == "@"):
                    dict_gre[ft] = {}
                    dict_gre[ft]["TF_NAME"] = data[1].rstrip()
                    dict_gre[ft]["SOURCE_NAME"] = data[2].rstrip()
                    dict_gre[ft]["TFBS_LENGHT"] = int(data[3].rstrip())
                    dict_gre[ft]["SAMPLE_SIZE"] = int(data[4].rstrip())
                    dict_gre[ft]["FT_CODE"] = data[5].rstrip()
                    dict_gre[ft]["SPECIE_NAME"] = data[6].rstrip()
                    continue
                elif (data[0] == "#"):
                    continue
                elif (data[0] == "*"):
                    dict_gre[ft][data[1]] = {}
                    dict_gre[ft][data[1]]["P1"] = float(data[2].rstrip())
                    dict_gre[ft][data[1]]["P2"] = float(data[3].rstrip())
                    dict_gre[ft][data[1]]["P3"] = float(data[4].rstrip())
                    dict_gre[ft][data[1]]["AUC mean"] = data[5].rstrip()
                    dict_gre[ft][data[1]]["Cutoff_0_1"] = float(data[6].rstrip())
                    dict_gre[ft][data[1]]["Cutoff_Youden"] = float(data[7].rstrip())
                    gre_file_name = tf_name + "_" + data[1].rstrip() + "_LAPFA_" + str(int(float(data[2].rstrip()))) + \
                                    "_" + str(data[3].rstrip()) + "_" + str(data[4].rstrip()) + ".bnf"
                    dict_gre[ft][data[1]]["GRE_FILENAME"] = gre_file_name
                    dict_gre[ft][data[1]]["GRE"] = {}
                    dict_gre[ft][data[1]]["GRE"] = parse_bnf(gre_file_name)
                    if (data[1] == "N3"):
                        break
    return dict_gre


def calc_logP_N1(seq):
    return len(seq) * math.log(0.25, 2)


def calc_logP_N2(seq, freq):
    result = 0.0
    for base in seq:
        if base.upper() == "A":
            result += math.log(freq['A'], 2)
        elif base.upper() == "C":
            result += math.log(freq["C"], 2)
        elif base.upper() == "G":
            result += math.log(freq["G"], 2)
        elif base.upper() == "T":
            result += math.log(freq["T"], 2)
    return result


def calc_logP_N3(seq):
    psc = global_pseudocontador_seq # 0.01
    result = 0.0
    for base in seq:
        if base == "A":
            result += math.log((seq.count('A') + psc) / len(seq), 2)
        elif base == "C":
            result += math.log((seq.count('C') + psc) / len(seq), 2)
        elif base == "G":
            result += math.log((seq.count('G') + psc) / len(seq), 2)
        elif base == "T":
            result += math.log((seq.count('T') + psc) / len(seq), 2)
    return result


def calc_score_lapfa(sliding_window, gre_file):

    os.chdir(global_compila_main)

    p = subprocess.Popen(global_compila_main + "/main " + gre_file,
                         shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line,
    retval = p.wait()

    p = subprocess.Popen("g++ Main2.cpp -o main2", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line,
    retval = p.wait()

    result = 0.0
    p = subprocess.Popen(global_compila_main + "/main2 " + sliding_window, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        print line
        result = float(line)
    retval = p.wait()

    return result


def parse_exe_par(par):
    dict_exe = {}
    with open(par, 'r') as exe_par_hdl:
        for line in exe_par_hdl:
            line = line.rstrip()
            if(line[0:1]=="#"):
                continue
            sep_pos = line.find("=")
            if (sep_pos == -1):
                continue
            left_line = line[0:sep_pos]
            right_line = line[sep_pos + 1:]
            words = right_line.split(",")
            try:
                dict_exe[left_line] = [int(i) for i in words]
            except ValueError:
                dict_exe[left_line] = words
    return dict_exe


def calc_score_pwm(sequence, pwm):
    score = 0.0
    i = 0
    for base in sequence:
        if base not in ('A', 'C', 'G', 'T'):
            continue
        score = score + pwm[base][i]
        i += 1
    return score


def calc_score_gre(tf, sequence, gre):
    knot = 1
    score = 0.0
    for i, base in enumerate(sequence):
        if base not in ["A","C","G","T"]:
            real_base = "#"
        else:
            real_base = base
        try:
            score += gre[knot][real_base][1]
            knot = gre[knot][real_base][0]
        except:
            global_gre_error.add(tf)
    score += gre[knot]["#"][1]
    return score


def get_number(str_number):
    first_pos = str_number.find("<")
    last_pos =  str_number.find(">")
    num_str = str_number[first_pos+1:last_pos]
    ret_number = int(num_str)
    return ret_number


def get_probability(str_number):
    first_pos = str_number.find("[")
    last_pos =  str_number.find("]")
    num_str = str_number[first_pos+1:last_pos]
    return float(num_str)


def parse_bnf(gre):
    os.chdir(global_gre_pwm_source)
    neper_number = exp(1)
    dict_gre = {}
    knot_son = 0
    with open(gre,'r') as gre_f:
        for line_gre in gre_f:
            sep_pos = line_gre.find("::=")
            if (sep_pos > 0):
                left_line  = line_gre[0:sep_pos]
                right_line = line_gre[sep_pos+3:]
                base_letter = right_line[1:2]
                knot_father = get_number(left_line)
                if base_letter != "#":
                    knot_son = get_number(right_line)
                if knot_father not in dict_gre:
                    dict_gre[knot_father] = {}
                dict_gre[knot_father][base_letter] = [knot_son, log(get_probability(right_line), neper_number)]
    return dict_gre


def ini_dict_count_tfbs(dict_exe_par):
    dict_count = {}
    for tf_id in dict_exe_par["FT_LIST"]:
        dict_count[tf_id] = {}
        dict_count[tf_id]["PWM"] = {}
        dict_count[tf_id]["GRE"] = {}
        for model in dict_exe_par["MODELS"]:
            dict_count[tf_id]["PWM"][model] = {}
            dict_count[tf_id]["GRE"][model] = {}
            for criterion in dict_exe_par["CRITERIA"]:
                dict_count[tf_id]["PWM"][model][criterion] = 0
                dict_count[tf_id]["GRE"][model][criterion] = 0
    return dict_count


def proc_fasta(exe_par):

    global global_compila_main, global_gre_pwm_source, global_pseudocontador_seq, global_gre_error
    global_gre_error = set()
    dict_exe_param = parse_exe_par(exe_par)
    global_compila_main = dict_exe_param["AUTOMATON_DIR"][0]
    global_gre_pwm_source = dict_exe_param["PWM_GRE_DIRECTORY_IN"][0]

    global_pseudocontador_seq = float(dict_exe_param["PSEUDOCOUNTER"][0])

    dict_pwm_parsed = parse_pwm(dict_exe_param["FT_LIST"], dict_exe_param["MODELS"])
    dict_gre_parsed = parse_gre(dict_exe_param["FT_LIST"], dict_exe_param["GRE_TOOL_PARAM"][0])
    dict_tf_param = parse_tf_par(dict_exe_param["FT_LIST"], dict_exe_param["TF_PARAM"][0])
    dict_freq_bas_especie = parse_freq_gen(dict_exe_param["GENOME_FREQ"][0])

    count_pos_block = -1
    count_pos_global = -1
    count_block = 1
    block_size = dict_exe_param["BLOCK_SIZE"][0]

    dict_count_tfbs = {}

    if (dict_exe_param["DEBUG"][0]=="NO"):
        console_print = False
    else:
        console_print = True

    fasta_sequences = SeqIO.parse(open(dict_exe_param["INPUT_DATA"][0]),'fasta')

    with open(dict_exe_param["SUMMARY_OUTPUT_DATA"][0], 'w') as summ_out_file:
        summ_out_file.write("METHOD"  + "\t" + "CROMOSSOME" + "\t" + "BLOCK#" + "\t" + "BLK_POS_INI" + "\t" +
                            "BLK_POS_END" + "\t" + "MODEL" + "\t" + "CRITERION" )
        for ddtf_id in dict_exe_param["FT_LIST"]:
            summ_out_file.write("\t" + "FT_" + str(ddtf_id))
        summ_out_file.write("\n")

        with open(dict_exe_param["OUTPUT_DATA"][0], 'w') as out_file:
            last_sliding_window = ""
            last_name = ""
            if (console_print):
                print("METHOD" + "\t" + "TF_CODE" + "\t" + "CROMOSSOME" + "\t"
                           + "GLOBAL_POSITION" + "\t" + "LOCAL_POSITION" + "\t" + "SLIDING_WINDOW" + "\t"
                           + "MODEL" + "\t" + "SCORE_PWM/GRE" + "\t" + "CUTOFF_PWM/GRE" + "\t" + "CRITERION" + "\t"
                           + "TF_ACRONYM" + "\t" + "SPECIE_NAME")
            out_file.write("METHOD" + "\t" + "TF_CODE" + "\t" + "CROMOSSOME" + "\t"
                           + "GLOBAL_POSITION" + "\t" + "LOCAL_POSITION" + "\t" + "SLIDING_WINDOW" + "\t"
                           + "MODEL" + "\t" + "SCORE_PWM/GRE" + "\t" + "CUTOFF_PWM/GRE" + "\t" + "CRITERION" + "\t"
                           + "TF_ACRONYM" + "\t" + "SPECIE_NAME" + "\n")

            dict_count_tfbs = ini_dict_count_tfbs(dict_exe_param)

            for fasta in fasta_sequences:
                sliding_window = ""
                cromossome, sequence = fasta.id, str(fasta.seq)
                if (cromossome not in dict_exe_param["CROMOSSOMES"]):
                    continue
                if (last_name == cromossome):
                    sequence = last_sliding_window + sequence
                len_seq = len(sequence)
                for i, base in enumerate(sequence):
                    count_pos_global += 1
                    count_pos_block += 1
                    for tf_id in dict_exe_param["FT_LIST"]:
                        tfbs_len = dict_tf_param["FT" + str(tf_id).zfill(7)]["TF_TFBS_LEN"]
                        sliding_window = sequence[i:i + tfbs_len]
                        if (i > len_seq-tfbs_len):
                            break
                        for model in dict_exe_param["MODELS"]:
                            for criterion in dict_exe_param["CRITERIA"]:
                                gre_cutoff = dict_gre_parsed[tf_id][model]["Cutoff_"+criterion]
                                pwm_cutoff = dict_pwm_parsed[tf_id][model]["Cutoff: "+criterion+" mean"]

                                score_pwm = calc_score_pwm(sliding_window, dict_pwm_parsed[tf_id][model]["PWM"])
                                score_gre = calc_score_gre(tf_id, sliding_window, dict_gre_parsed[tf_id][model]["GRE"])

                                if (console_print):
                                    score_lapfa = calc_score_lapfa(sliding_window, global_gre_pwm_source +
                                                dict_gre_parsed[tf_id][model]["GRE_FILENAME"])
                                    print("Escore automato GrammarLab =" + str(score_lapfa)  +
                                          "  | Escore calculado =" + str(score_gre) +
                                          "   (antes de descontar a freq. de fundo).")

                                if (model == "N1"):
                                    score_gre = score_gre - calc_logP_N1(sliding_window)
                                elif (model == "N2"):
                                    score_gre = score_gre - calc_logP_N2(sliding_window,
                                                        dict_freq_bas_especie[dict_gre_parsed[tf_id]["SPECIE_NAME"]])
                                elif (model == "N3"):
                                    score_gre = score_gre - calc_logP_N3(sliding_window)

                                if score_pwm > pwm_cutoff:
                                    out_file.write("PWM" + "\t" + "FT" + str(tf_id).zfill(7) + "\t" + cromossome + "\t"
                                    + str(count_pos_global) + "\t" + str(i) + "\t" + sliding_window + "\t"
                                    + model + "\t" + str(score_pwm) + "\t" + str(pwm_cutoff) + "\t" + criterion + "\t"
                                    + dict_gre_parsed[tf_id]["FT_CODE"] + "\t" + dict_gre_parsed[tf_id]["SPECIE_NAME"]
                                    + "\n")

                                    dict_count_tfbs[tf_id]["PWM"][model][criterion] += 1

                                if score_gre > gre_cutoff:
                                    out_file.write("GRE" + "\t" + "FT" + str(tf_id).zfill(7) + "\t" + cromossome + "\t"
                                    + str(count_pos_global) + "\t" + str(i) + "\t" + sliding_window + "\t"
                                    + model + "\t" + str(score_gre) + "\t" + str(gre_cutoff) + "\t" + criterion + "\t"
                                    + dict_gre_parsed[tf_id]["FT_CODE"] + "\t" + dict_gre_parsed[tf_id]["SPECIE_NAME"]
                                    + "\n")

                                    dict_count_tfbs[tf_id]["GRE"][model][criterion] += 1

                                if (count_pos_block >= block_size):
                                    block_range = str(count_pos_global-block_size) + "\t" + str(count_pos_global)

                                    for md_id in dict_exe_param["MODELS"]:
                                        for cr_id in dict_exe_param["CRITERIA"]:
                                            summ_out_file.write("PWM:" + "\t" + cromossome + "\t" + str(count_block) +
                                            "\t" + block_range + "\t" + md_id + "\t"  + cr_id )
                                            for mtf_id in dict_exe_param["FT_LIST"]:
                                                summ_out_file.write("\t" + str(dict_count_tfbs[mtf_id]["PWM"][md_id][cr_id]))
                                            summ_out_file.write("\n")

                                    for md_id in dict_exe_param["MODELS"]:
                                        for cr_id in dict_exe_param["CRITERIA"]:
                                            summ_out_file.write("GRE:" + "\t" + cromossome + "\t" + str(count_block) +
                                                                "\t" + block_range + "\t" + md_id + "\t" + cr_id)
                                            for mtf_id in dict_exe_param["FT_LIST"]:
                                                summ_out_file.write(
                                                    "\t" + str(dict_count_tfbs[mtf_id]["GRE"][md_id][cr_id]))
                                            summ_out_file.write("\n")

                                    dict_count_tfbs = ini_dict_count_tfbs(dict_exe_param)
                                    count_block += 1
                                    count_pos_block = -1

                                    out_file.flush()
                                    summ_out_file.flush()

                last_sliding_window = sliding_window
                last_name = cromossome

    print ("GREs found with errors:")
    print " ".join(str(x) for x in global_gre_error)


proc_fasta("/home/aferrao/tfbs_search_tool/exe_parameters.txt")

