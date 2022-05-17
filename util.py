import os

def mkdirs(path):
    """
    Generate folder (including the path) based on the input path

    Parameters
    ----------
    path : str
        The input path to be created
    """

    isExists = os.path.exists(path)

    if not isExists:
        os.makedirs(path)

def merge_file(input_fn_list, output_fn):
    """concatenate files in the file list to generate a new single file.

    Parameters
    ----------
    input_fn_list : list
        The list of files to be concatenated
    output_fn : str
        The output file path + name 
    """
    file_list_temp = " ".join(input_fn_list)
    print("Merging files:", file_list_temp)
    os.system("cat " + file_list_temp + " > " + output_fn)

if __name__ == '__main__':
    # path = "./result/uniprot-id_25/policy2/"
    # dataset = "uniprot-id_25"

    path = "./result/uniprot-human_proteome-reviewed/policy2/"
    dataset = "uniprot-human_proteome-reviewed"

    # file_list = [path + "uniprot-id_25_output_C_KR.fasta", path + "uniprot-id_25_output_E_KR.fasta"]

    # merge_file(file_list, path + "uniprot-id_25_output_C+E_KR.fasta")

    file_list = [path + dataset +"_output_AllExceptC_KR.fasta", path + dataset + "_output_C2B_KR.fasta", path + dataset + "_output_KR_KR_C2Z.fasta", path + dataset + "_output_KR_KR_D2U.fasta", path + dataset + "_output_KR_KR_E2O.fasta"]

    merge_file(file_list, path + dataset + "_output_AllExceptC+C2B_trans_KR.fasta")