from Bio import SeqIO
from util import *
from Bio.SeqRecord import SeqRecord

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
        elif elem in end_elems:
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



if __name__ == '__main__':
    input_dict = SeqIO.index_db("./dataset/uniprot-id_25/uniprot-id_25.idx", "./dataset/uniprot-id_25/uniprot-id_25.fasta", "fasta")

    record = input_dict["sp|O75452|RDH16_HUMAN"]
    print("record.seq:")
    print(record.seq)
    print("record.id", record.id)
    print("record.description", record.description)

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
    print("new_description", new_description)