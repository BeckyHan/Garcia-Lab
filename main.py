from fasta_proc_lib import *
import copy


if __name__ == '__main__':
    input_fn = "./dataset/uniprot-id_25/uniprot-id_25.fasta"
    
    # out_path = "./result"
    # out_fn = "uniprot-id_25_output.fasta"

    # input_fn_idx = "./dataset/uniprot-human_proteome-reviewed/uniprot-human_proteome-reviewed.idx"
    # input_fn = "./dataset/uniprot-human_proteome-reviewed/uniprot-human_proteome-reviewed.fasta"
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
    print("input_path:", input_path)
    print("input_fn_idx", input_fn_idx)
    print("output_path:", output_path)
    
    
    input_dict = SeqIO.index_db(input_fn_idx, input_fn, "fasta")

    rec_output_all = []
    for key in input_dict.keys():
        rec_output_lst = fasta_proc_1(None, input_dict[key])
        rec_output_all += rec_output_lst
    
    mkdirs(output_path)
    SeqIO.write(rec_output_all, output_path + "/" + file_name + "_D_output.fasta", "fasta")

