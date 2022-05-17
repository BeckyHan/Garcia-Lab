from tkinter.messagebox import NO
from Bio import SeqIO
from util import *
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def fasta_proc_1(opts, record, start_elems = ["D"], end_elems = ["K", "R"]):
    """extract a subsequence from the sequence of seq_record based on the rule and save them in files.
    
    Starts from start_elems and end at end_elems. By default, start_elems = ["D"] and end_elems = ["K", "R"]. Here we use this default case as an example

    1. Subseq starts from D and ends at K or R.
    2. Insert R in front of the beginnning D
    3. Extract multiple subsubseq from the subseq based on the same rule if there are multiple D in the subseq
    4. If there is R in front of D in the original seq, the save file name with _R
    5. Extract multiple subsubseq from the subseq based on the same rule if there are multiple K/R at the end of subseq (e.g., XXXKRK)


    Parameters
    ----------
    opts : class
        input options class
    record : Bio.SeqRecord.SeqRecord
        a Biopy SeqRecord to be processed
    start_elems : list
        list of start elements in the subsequence. "D" by default
    end_elems : list
        list of end elements in the subsequence. "K" or "R" by default
    """
    # genrate a part of subseq ID
    start_elem_name = "".join(start_elems)
    end_elem_name = "".join(end_elems)

    # generate the new description for the newly generated subsequence
    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    new_description = " ".join(new_description_temp)

    seq = record.seq
    # find index of D, K and R in original Seq
    D_idx_lst = []
    KR_idx_lst = []

    for elem_idx, elem in enumerate(seq):
        if elem in start_elems:
            D_idx_lst.append(elem_idx)
        if elem in end_elems:
            KR_idx_lst.append(elem_idx)
    
    subseq_idx = 0 # used to generate the subseq file name
    rec_output_lst = []
    if len(D_idx_lst) != 0 and len(KR_idx_lst) != 0:
        for D_idx in D_idx_lst:
            # make sure there is a KR_idx so that KR_idx > D_idx
            if (KR_idx_lst[-1] < D_idx):
                break
            else:
                # find the first KR_idx_temp so that KR_idx_temp > D_idx
                for KR_idx in KR_idx_lst:
                    if KR_idx > D_idx:
                        KR_idx_temp = KR_idx
                        break
                
                # check and generate subseq alone with seq[KR_idx_temp], seq[KR_idx_temp + 1], ..., seq[len(seq) - 1]
                # process the case in which K/R is all together (e.g., XXXKRK)
                for i in range(KR_idx_temp, len(seq)):
                    if seq[i] in end_elems:
                        subseq = seq[D_idx : i] + seq[i]
                        # Add R in the front of subseq
                        subseq = "R" + subseq

                        # generate the file_name for the subseq
                        subseq_name = "subseq%d" % subseq_idx + "_" + start_elem_name + "_" + end_elem_name
                        if D_idx != 0 and seq[D_idx - 1] == "R":
                            subseq_name = subseq_name + "_R"
                        subseq_idx += 1
                        
                        # with open(full_path + "/" + subseq_name + "txt", 'w') as f:
                        #     f.write(str(subseq))

                        rec_output_lst.append(SeqRecord(subseq, id=record.id+"_"+subseq_name, description=new_description))

                    else:
                        break
    return rec_output_lst

def fasta_proc_2(opts, record, start_elems = ["D"], end_elems = ["K", "R"]):
    """extract a subsequence from the sequence of seq_record based on the rule and save them in files.
    
    Starts from start_elems and end at end_elems. By default, start_elems = ["D"] and end_elems = ["K", "R"]. Here we use this default case as an example

    1. Subseq starts from D and ends at K or R.
    2. Insert R in front of the beginnning D
    3. Extract multiple subsubseq from the subseq based on the same rule if there are multiple D in the subseq
    4. If there is R in front of D in the original seq, the save file name with _R
    5. Extract multiple subsubseq from the subseq based on the same rule if there are multiple K/R at the end of subseq (e.g., XXXKRK)
    6. ignore subseqs whose length <= 6
    7. Remove a part of description of all the sebseqs behind "OS"
    8. Besides above rules, a different rule (miss damage) is also added that DXXXK/RXXXK/R is also considered

    Parameters
    ----------
    opts : class
        input options class
    record : Bio.SeqRecord.SeqRecord
        a Biopy SeqRecord to be processed
    start_elems : list
        list of start elements in the subsequence. "D" by default
    end_elems : list
        list of end elements in the subsequence. "K" or "R" by default
    
    Returns
    ----------
    rec_output_lst : list
        list of Subseq Record
    """
    # genrate a part of subseq ID
    start_elem_name = "".join(start_elems)
    end_elem_name = "".join(end_elems)

    # generate the new description for the newly generated subsequence
    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    new_description = " ".join(new_description_temp)
    new_description_temp = new_description.split("OS")
    new_description_temp.pop(-1)
    new_description = " ".join(new_description_temp)

    seq = record.seq
    # find index of D, K and R in original Seq
    D_idx_lst = []
    KR_idx_lst = []

    for elem_idx, elem in enumerate(seq):
        if elem in start_elems:
            D_idx_lst.append(elem_idx)
        if elem in end_elems:
            KR_idx_lst.append(elem_idx)
    
    subseq_idx = 0 # used to generate the subseq file name
    rec_output_lst = []
    if len(D_idx_lst) != 0 and len(KR_idx_lst) != 0:

        # Policy 1
        for D_idx in D_idx_lst:
            # make sure there is a KR_idx so that KR_idx > D_idx
            if (KR_idx_lst[-1] <= D_idx):
                break
            else:
                # find the first KR_idx_temp so that KR_idx_temp > D_idx
                for KR_idx in KR_idx_lst:
                    if KR_idx > D_idx:
                        KR_idx_temp = KR_idx
                        break
                
                # check and generate subseq alone with seq[KR_idx_temp], seq[KR_idx_temp + 1], ..., seq[len(seq) - 1]
                # process the case in which K/R is all together (e.g., XXXKRK)
                for i in range(KR_idx_temp, len(seq)):
                    if seq[i] in end_elems:
                        subseq = seq[D_idx : i] + seq[i]
                        # Add R in the front of subseq
                        subseq = "R" + subseq

                        # ignore the subseq whose length is less than 6
                        if len(subseq) <= 6:
                            continue

                        # generate the file_name for the subseq
                        subseq_name = "subseq%d" % subseq_idx + "_" + start_elem_name + "_" + end_elem_name
                        if D_idx != 0 and seq[D_idx - 1] == "R":
                            subseq_name = subseq_name + "_R"
                        subseq_idx += 1
                        
                        # with open(full_path + "/" + subseq_name + "txt", 'w') as f:
                        #     f.write(str(subseq))

                        rec_output_lst.append(SeqRecord(subseq, id=record.id+"_"+subseq_name, description=new_description))

                    else:
                        break
        
        # Policy 2 (miss damage)
        for D_idx in D_idx_lst:
            # make sure there is a KR_idx so that KR_idx > D_idx
            if (KR_idx_lst[-1] <= D_idx):
                break
            else:
                # # find the first KR_idx_temp so that KR_idx_temp > D_idx
                # for KR_idx in KR_idx_lst:
                #     if KR_idx > D_idx:
                #         KR_idx_temp = KR_idx
                #         break

                # find the first KR_idx_temp so that KR_idx_temp > D_idx and KR_idx_temp is either the first KK/KR/RK/RR... or the second K/R
                cnt = 0
                for KR_idx in KR_idx_lst:
                    if KR_idx > D_idx:
                        # KR_idx_temp is the first KK/KR/RK/RR...
                        # Since this case has been processed in policy 1, ignore it in this policy
                        if (KR_idx + 1) in KR_idx_lst:
                                # KR_idx_temp = KR_idx
                                KR_idx_temp = len(seq)
                                break
                        # or KR_idx_temp is the second K/R
                        elif cnt == 1:
                            KR_idx_temp = KR_idx
                            break
                        else:
                            pass
                        
                        # terminate the whole program if subseq D....K/R in which K/R is the last one in seq
                        if KR_idx == KR_idx_lst[-1] and cnt == 0:
                            KR_idx_temp = len(seq)
                            break
                        
                        cnt += 1
                        
                # check and generate subseq alone with seq[KR_idx_temp], seq[KR_idx_temp + 1], ..., seq[len(seq) - 1]
                # process the case in which K/R is all together (e.g., XXXKRK)
                for i in range(KR_idx_temp, len(seq)):
                    if seq[i] in end_elems:
                        subseq = seq[D_idx : i] + seq[i]
                        # Add R in the front of subseq
                        subseq = "R" + subseq

                        # ignore the subseq whose length is less than 6
                        if len(subseq) <= 6:
                            continue

                        # generate the file_name for the subseq
                        subseq_name = "subseq%d" % subseq_idx + "_" + start_elem_name + "_" + end_elem_name
                        if D_idx != 0 and seq[D_idx - 1] == "R":
                            subseq_name = subseq_name + "_R"
                        subseq_idx += 1
                        
                        # with open(full_path + "/" + subseq_name + "txt", 'w') as f:
                        #     f.write(str(subseq))

                        rec_output_lst.append(SeqRecord(subseq, id=record.id+"_"+subseq_name, description=new_description))

                    else:
                        break
    return rec_output_lst

def fasta_proc_3(opts, record, start_elems = ["D"], end_elems = ["K", "R"]):
    """extract a subsequence from the sequence of seq_record based on the rule and save them in files.
    
    Starts from start_elems and end at end_elems. By default, start_elems = ["D"] and end_elems = ["K", "R"]. Here we use this default case as an example. The extracted subseq should be D|XXXXK/R (starting elem D is not included in subseq).

    1. Subseq starts from D and ends at K or R.
    2. Remove start elems (D by default) from the subseq
    3. Extract multiple subsubseq from the subseq based on the same rule if there are multiple start elems in the subseq
    4. If there is R in front of D in the original seq, the save file name with _R
    5. Extract multiple subsubseq from the subseq based on the same rule if there are multiple K/R at the end of subseq (e.g., XXXKRK)
    6. ignore subseqs whose length <= 5
    7. Remove a part of description of all the sebseqs behind "OS"

    Parameters
    ----------
    opts : class
        input options class
    record : Bio.SeqRecord.SeqRecord
        a Biopy SeqRecord to be processed
    start_elems : list
        list of start elements in the subsequence. "D" by default
    end_elems : list
        list of end elements in the subsequence. "K" or "R" by default
    
    Returns
    ----------
    rec_output_lst : list
        list of Subseq Record
    """
    # genrate a part of subseq ID
    start_elem_name = "".join(start_elems)
    end_elem_name = "".join(end_elems)

    # generate the new description for the newly generated subsequence
    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    new_description = " ".join(new_description_temp)
    new_description_temp = new_description.split("OS")
    new_description_temp.pop(-1)
    new_description = " ".join(new_description_temp)

    seq = record.seq
    # find index of D, K and R in original Seq
    D_idx_lst = []
    KR_idx_lst = []

    for elem_idx, elem in enumerate(seq):
        if elem in start_elems:
            D_idx_lst.append(elem_idx)
        if elem in end_elems:
            KR_idx_lst.append(elem_idx)
    
    subseq_idx = 0 # used to generate the subseq file name
    rec_output_lst = []
    if len(D_idx_lst) != 0 and len(KR_idx_lst) != 0:

        # Policy 1
        for D_idx in D_idx_lst:
            # make sure there is a KR_idx so that KR_idx > D_idx
            if (KR_idx_lst[-1] <= D_idx):
                break
            else:
                # find the first KR_idx_temp so that KR_idx_temp > D_idx
                for KR_idx in KR_idx_lst:
                    if KR_idx > D_idx:
                        KR_idx_temp = KR_idx
                        break
                
                # check and generate subseq alone with seq[KR_idx_temp], seq[KR_idx_temp + 1], ..., seq[len(seq) - 1]
                # process the case in which K/R is all together (e.g., XXXKRK)
                for i in range(KR_idx_temp, len(seq)):
                    if seq[i] in end_elems:
                        subseq = seq[D_idx+1 : i] + seq[i]
                        # Add R in the front of subseq
                        # subseq = "R" + subseq

                        # ignore the subseq whose length is less than 6
                        if len(subseq) <= 6:
                            continue

                        # generate the file_name for the subseq
                        subseq_name = "subseq%d" % subseq_idx + "_" + start_elem_name + "_" + end_elem_name
                        if D_idx != 0 and seq[D_idx - 1] == "R":
                            subseq_name = subseq_name + "_R"
                        subseq_idx += 1
                        
                        # with open(full_path + "/" + subseq_name + "txt", 'w') as f:
                        #     f.write(str(subseq))

                        rec_output_lst.append(SeqRecord(subseq, id=record.id+"_"+subseq_name, description=new_description))

                    else:
                        break
        
        # Policy 2 (miss damage)
        for D_idx in D_idx_lst:
            # make sure there is a KR_idx so that KR_idx > D_idx
            if (KR_idx_lst[-1] <= D_idx):
                break
            else:
                # # find the first KR_idx_temp so that KR_idx_temp > D_idx
                # for KR_idx in KR_idx_lst:
                #     if KR_idx > D_idx:
                #         KR_idx_temp = KR_idx
                #         break

                # find the first KR_idx_temp so that KR_idx_temp > D_idx and KR_idx_temp is either the first KK/KR/RK/RR... or the second K/R
                cnt = 0
                for KR_idx in KR_idx_lst:
                    if KR_idx > D_idx:
                        # KR_idx_temp is the first KK/KR/RK/RR...
                        # Since this case has been processed in policy 1, ignore it in this policy
                        if (KR_idx + 1) in KR_idx_lst:
                                # KR_idx_temp = KR_idx
                                KR_idx_temp = len(seq)
                                break
                        # or KR_idx_temp is the second K/R
                        elif cnt == 1:
                            KR_idx_temp = KR_idx
                            break
                        else:
                            pass
                        
                        # terminate the whole program if subseq D....K/R in which K/R is the last one in seq
                        if KR_idx == KR_idx_lst[-1] and cnt == 0:
                            KR_idx_temp = len(seq)
                            break
                        
                        cnt += 1
                        
                # check and generate subseq alone with seq[KR_idx_temp], seq[KR_idx_temp + 1], ..., seq[len(seq) - 1]
                # process the case in which K/R is all together (e.g., XXXKRK)
                for i in range(KR_idx_temp, len(seq)):
                    if seq[i] in end_elems:
                        subseq = seq[D_idx+1 : i] + seq[i]
                        # Add R in the front of subseq
                        # subseq = "R" + subseq

                        # ignore the subseq whose length is less than 6
                        if len(subseq) <= 6:
                            continue

                        # generate the file_name for the subseq
                        subseq_name = "subseq%d" % subseq_idx + "_" + start_elem_name + "_" + end_elem_name
                        if D_idx != 0 and seq[D_idx - 1] == "R":
                            subseq_name = subseq_name + "_R"
                        subseq_idx += 1
                        
                        # with open(full_path + "/" + subseq_name + "txt", 'w') as f:
                        #     f.write(str(subseq))

                        rec_output_lst.append(SeqRecord(subseq, id=record.id+"_"+subseq_name, description=new_description))

                    else:
                        break
    return rec_output_lst


def fasta_trans_first_elem(opts, record, start_elems = ["C"], replace_elem = "B"):
    """replace the first C in the input seq starting with "RC..." with "B", so that the seq starts with "RB...."

    Parameters
    ----------
    opts : class
        input options class
    record : Bio.SeqRecord.SeqRecord
        a Biopy SeqRecord to be processed
    start_elems : list
        list of start elements in the seq. "C" by default
    replace_elem : str
        Elem used to replace the first start element in seq. "B" by default

    Returns
    ----------
    record_new : Record
        new record with modified seq
    """
    # generate the new description for the newly generated subsequence
    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    new_description = " ".join(new_description_temp)

    # generate new id for each seq
    # replace the starting element to replace policy
    # e.g. replace "XXX_C_XXX" to "XXX_C2B_XXX"
    new_id_elem = "".join(start_elems) + "2" + replace_elem
    id_list_temp = record.id.split("_")
    id_list_temp2 = []
    for id_elem in id_list_temp:
        if id_elem in start_elems:
            id_elem = new_id_elem
        id_list_temp2.append(id_elem)
    new_id = "_".join(id_list_temp2)

    seq = record.seq

    target_elem = seq[1]
    if target_elem not in start_elems:
        raise ValueError("target element is not contained in start_elems!", target_elem, start_elems)
    
    seq_new_lst = list(seq)
    seq_new_lst[1] = replace_elem
    seq_new = "".join(seq_new_lst)
    record_new = SeqRecord(Seq(seq_new), id=new_id, description=new_description)

    return record_new

def fasta_trans_elem_once(opts, record, trans_dict):
    """Find all the target elems in the seq of record and replace it with replace_elem. One newly generated record/seq for each translation

    Parameters
    ----------
    opts : class
        input options class
    record : Bio.SeqRecord.SeqRecord
        a Biopy SeqRecord to be processed
    trans_dict : dictionary
        Translation dictionary in which key represents the targets elements to be replaced and vaule represents the corresponding replaced elements (target_elem : replace_elem) 

    Returns
    ----------
    rec_output_lst : list
        list of Subseq Record
    """
    # generate the new description for the newly generated subsequence
    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    new_description = " ".join(new_description_temp)

    # generate new id segment for each seq
    id_seg_list_temp = [key + "2" + trans_dict[key] for key in trans_dict.keys()]
    new_id_seg = "_".join(id_seg_list_temp)


    seq = record.seq

    subseq_idx = 0 # used to generate the subseq file name
    rec_output_lst = []

    for i in range(0, len(seq)):
        replace_elem = trans_dict.get(seq[i])
        if replace_elem:
            seq_new_lst = list(seq)
            seq_new_lst[i] = replace_elem
            seq_new = "".join(seq_new_lst)
            record_new = SeqRecord(Seq(seq_new), id=record.id+"_"+new_id_seg+"_"+str(subseq_idx), description=new_description)
            rec_output_lst.append(record_new)
            subseq_idx += 1

    return rec_output_lst




if __name__ == '__main__':
    input_dict = SeqIO.index_db("./dataset/uniprot-id_25/uniprot-id_25.idx", "./dataset/uniprot-id_25/uniprot-id_25.fasta", "fasta")

    record = input_dict["sp|O75452|RDH16_HUMAN"]
    print("record.seq:")
    print(record.seq)
    print("record.id", record.id)
    print("record.description", record.description)

    # **************** Policy1 ****************
    # calculate for subseq starting from D and end at K or R
    rec_output_lst = fasta_proc_1(None, record, start_elems = ["D"], end_elems = ["K", "R"])
    print(len(rec_output_lst))

    mkdirs("./test")
    SeqIO.write(rec_output_lst,"./test/uniprot-id_25_output_D_KR.fasta", "fasta")

    # calculate for subseq starting from E and end at K or R
    rec_output_lst = fasta_proc_1(None, record, start_elems = ["E"], end_elems = ["K", "R"])
    print(len(rec_output_lst))

    mkdirs("./test")
    SeqIO.write(rec_output_lst,"./test/uniprot-id_25_output_E_KR.fasta", "fasta")

    # calculate for subseq starting from C and end at K or R
    rec_output_lst = fasta_proc_1(None, record, start_elems = ["C"], end_elems = ["K", "R"])
    print(len(rec_output_lst))

    mkdirs("./test")
    SeqIO.write(rec_output_lst,"./test/uniprot-id_25_output_C_KR.fasta", "fasta")

    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    print(new_description_temp)
    new_description = " ".join(new_description_temp)
    print("new_description:", new_description)

    # **************** Policy2 ****************
    print("Performing Policy2...")
    # calculate for subseq starting from D and end at K or R
    rec_output_lst = fasta_proc_2(None, record, start_elems = ["D"], end_elems = ["K", "R"])
    print(len(rec_output_lst))

    mkdirs("./test/policy2")
    SeqIO.write(rec_output_lst,"./test/policy2/uniprot-id_25_output_D_KR.fasta", "fasta")

    new_description_temp = record.description.split(" ")
    new_description_temp.pop(0)
    new_description = " ".join(new_description_temp)
    new_description_temp = new_description.split("OS")
    new_description_temp.pop(-1)
    print(new_description_temp)
    new_description = " ".join(new_description_temp)
    print("new_description:", new_description)

    # convert first C to B in XXX_output_C_KR.fasta
    print("Convert first C to B for all seq...")

    input_dict_c2b = SeqIO.index_db("./test/uniprot-id_25_output_C_KR.idx", "./test/uniprot-id_25_output_C_KR.fasta", "fasta")

    rec_convert_lst = []
    for key in input_dict_c2b.keys():
        rec_convert_lst_temp = fasta_trans_first_elem(None, input_dict_c2b[key], start_elems = ["C"], replace_elem = "B")
        rec_convert_lst.append(rec_convert_lst_temp)
    
    SeqIO.write(rec_convert_lst, "./test/uniprot-id_25_output_C2B_KR.fasta", "fasta")

    # calculate for subseq starting from K or R and end at K or R
    print("Converting C D E...")
    rec_output_lst = fasta_proc_3(None, record, start_elems = ["K", "R"], end_elems = ["K", "R"])
    SeqIO.write(rec_output_lst, "./test/policy2/uniprot-id_25_output_KR_KR.fasta", "fasta")

    # translate target elements in XXX_output_KR_KR.fasta. One record for each translation
    input_dict_trans = SeqIO.index_db("./test/policy2/uniprot-id_25_output_KR_KR.idx", "./test/policy2/uniprot-id_25_output_KR_KR.fasta", "fasta")

    record = input_dict_trans["sp|O75452|RDH16_HUMAN_subseq17_KR_KR"]
    print("record.seq:")
    print(record.seq)

    # trans_dict = {"C": "Z", "D": "U", "E": "O"}
    # rec_output_lst = fasta_trans_elem_once(None, record, trans_dict)
    # SeqIO.write(rec_output_lst, "./test/policy2/uniprot-id_25_output_KR_KR_trans.fasta", "fasta")

    trans_dict1 = {"C": "Z"}
    rec_output_lst = fasta_trans_elem_once(None, record, trans_dict1)
    SeqIO.write(rec_output_lst, "./test/policy2/uniprot-id_25_output_KR_KR_C2Z.fasta", "fasta")

    trans_dict2 = {"D": "U"}
    rec_output_lst = fasta_trans_elem_once(None, record, trans_dict2)
    SeqIO.write(rec_output_lst, "./test/policy2/uniprot-id_25_output_KR_KR_D2U.fasta", "fasta")

    trans_dict3 = {"E": "O"}
    rec_output_lst = fasta_trans_elem_once(None, record, trans_dict3)
    SeqIO.write(rec_output_lst, "./test/policy2/uniprot-id_25_output_KR_KR_E2O.fasta", "fasta")


