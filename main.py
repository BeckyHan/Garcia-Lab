from fasta_proc_lib import *
import copy


if __name__ == '__main__':
    # input_fn = "./dataset/uniprot-id_25/uniprot-id_25.fasta"
    input_fn = "./dataset/uniprot-human_proteome-reviewed/uniprot-human_proteome-reviewed.fasta"
    
    # out_path = "./result"
    # out_fn = "uniprot-id_25_output.fasta"
    
    # out_path = "./result"
    # out_fn = "uniprot-human_proteome-reviewed_output.fasta"

    # generate fn path
    fn_temp_lst1 = input_fn.split("/")
    file_name = fn_temp_lst1[-1].split(".").pop(0)
    fn_temp_lst1.pop(-1)
    input_path = "/".join(fn_temp_lst1)

    fn_temp_lst2 = copy.deepcopy(fn_temp_lst1)
    fn_temp_lst2[1] = "result"
    output_path = "/".join(fn_temp_lst2)

    # input_fn_idx
    input_fn_idx = input_path + "/" + file_name + ".idx"

    print("file_name:", file_name)
    print("input path:", input_path)
    print("input idx path", input_fn_idx)
    print("output path:", output_path)
    
    
    input_dict = SeqIO.index_db(input_fn_idx, input_fn, "fasta")

    # **************** Policy1 ****************
    # print("Performing Policy1...")
    # # calculate for subseq starting from D and end at K or R
    # rec_output_all = []
    # for key in input_dict.keys():
    #     rec_output_lst = fasta_proc_1(None, input_dict[key], start_elems = ["D"], end_elems = ["K", "R"])
    #     rec_output_all += rec_output_lst
    
    # mkdirs(output_path)
    # SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_D_KR.fasta", "fasta")

    # # calculate for subseq starting from E and end at K or R
    # rec_output_all = []
    # for key in input_dict.keys():
    #     rec_output_lst = fasta_proc_1(None, input_dict[key], start_elems = ["E"], end_elems = ["K", "R"])
    #     rec_output_all += rec_output_lst
    
    # mkdirs(output_path)
    # SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_E_KR.fasta", "fasta")

    # # calculate for subseq starting from C and end at K or R
    # rec_output_all = []
    # for key in input_dict.keys():
    #     rec_output_lst = fasta_proc_1(None, input_dict[key], start_elems = ["C"], end_elems = ["K", "R"])
    #     rec_output_all += rec_output_lst
    
    # mkdirs(output_path)
    # SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_C_KR.fasta", "fasta")

    # **************** Policy2 ****************
    print("Performing Policy2...")
    output_path = output_path + "/policy2"
    """
    # calculate for subseq starting from D and end at K or R
    rec_output_all = []
    for key in input_dict.keys():
        rec_output_lst = fasta_proc_2(None, input_dict[key], start_elems = ["D"], end_elems = ["K", "R"])
        rec_output_all += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_D_KR.fasta", "fasta")

    # calculate for subseq starting from E and end at K or R
    rec_output_all = []
    for key in input_dict.keys():
        rec_output_lst = fasta_proc_2(None, input_dict[key], start_elems = ["E"], end_elems = ["K", "R"])
        rec_output_all += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_E_KR.fasta", "fasta")

    # calculate for subseq starting from C and end at K or R
    rec_output_all = []
    for key in input_dict.keys():
        rec_output_lst = fasta_proc_2(None, input_dict[key], start_elems = ["C"], end_elems = ["K", "R"])
        rec_output_all += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_C_KR.fasta", "fasta")
    """
    # convert first C to B in XXX_output_C_KR.fasta
    input_dict_c2b = SeqIO.index_db(output_path + "/" + file_name + "_output_C_KR.idx", output_path + "/" + file_name + "_output_C_KR.fasta", "fasta")

    rec_convert_lst = []
    for key in input_dict_c2b.keys():
        rec_convert_temp = fasta_trans_first_elem(None, input_dict_c2b[key], start_elems = ["C"], replace_elem = "B")
        rec_convert_lst.append(rec_convert_temp)
    
    SeqIO.write(rec_convert_lst, output_path + "/" + file_name + "_output_C2B_KR.fasta", "fasta")


    # # calculate for subseq starting from K or R and end at K or R
    # rec_output_all = []
    # for key in input_dict.keys():
    #     rec_output_lst = fasta_proc_3(None, input_dict[key], start_elems = ["K", "R"], end_elems = ["K", "R"])
    #     rec_output_all += rec_output_lst
    
    # mkdirs(output_path)
    # SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_KR_KR.fasta", "fasta")

    # translate target elements in XXX_output_KR_KR.fasta
    # One record for each translation
    input_dict_trans = SeqIO.index_db(output_path + "/" + file_name + "_output_KR_KR.idx", output_path + "/" + file_name + "_output_KR_KR.fasta", "fasta")

    # translation table 1 - C to Z
    trans_dict1 = {"C": "Z"}

    rec_output_all_trans = []
    for key in input_dict_trans.keys():
        rec_output_lst = fasta_trans_elem_once(None, input_dict_trans[key], trans_dict1)
        rec_output_all_trans += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all_trans, output_path + "/" + file_name + "_output_KR_KR_C2Z.fasta", "fasta")

    # translation table 2 - D to U
    trans_dict2 = {"D": "U"}

    rec_output_all_trans = []
    for key in input_dict_trans.keys():
        rec_output_lst = fasta_trans_elem_once(None, input_dict_trans[key], trans_dict2)
        rec_output_all_trans += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all_trans, output_path + "/" + file_name + "_output_KR_KR_D2U.fasta", "fasta")

    # translation table 2 - E to O
    trans_dict3 = {"E": "O"}

    rec_output_all_trans = []
    for key in input_dict_trans.keys():
        rec_output_lst = fasta_trans_elem_once(None, input_dict_trans[key], trans_dict3)
        rec_output_all_trans += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all_trans, output_path + "/" + file_name + "_output_KR_KR_E2O.fasta", "fasta")


    # # calculate for subseq starting from all elements except for C and end at K or R
    # start_elem_list = ["A", "R", "N", "D", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    # rec_output_all = []
    # for key in input_dict.keys():
    #     rec_output_lst = fasta_proc_2(None, input_dict[key], start_elems = start_elem_list, end_elems = ["K", "R"])
    #     rec_output_all += rec_output_lst
    
    # mkdirs(output_path)
    # SeqIO.write(rec_output_all, output_path + "/" + file_name + "_output_AllExceptC_KR.fasta", "fasta")