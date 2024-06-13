import os
import sys
import concurrent.futures
from shutil import copyfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def create_fasta_files(contig_ids, lr_ids, contig_fasta, lr_fasta, output_dir, folder_num):
    folder_name = str(folder_num)
    output_dir = os.path.join(output_dir, folder_name)
    os.makedirs(output_dir, exist_ok=True)
    
    # Process contig sequences
    contig_out_file = os.path.join(output_dir, f"{folder_name}_Contigs.fasta")
    contig_records = SeqIO.to_dict(SeqIO.parse(contig_fasta, "fasta"))
    contig_seqs = []
    for contig_id in contig_ids:
        if contig_id in contig_records:
            contig_seq = contig_records[contig_id].seq
            contig_seqs.append(SeqRecord(contig_seq, id=contig_id, description=""))
        else:
            print(f"Warning: Contig {contig_id} not found in the contig FASTA file.")
    
    #with open(contig_out_file, 'w') as contig_fasta_file:
        #SeqIO.write(contig_seqs, contig_fasta_file, "fasta")
    
    # Copy contig_out_file to phase1_2_output.fasta
    phase1_2_output_file = os.path.join(input_dir, "phase1_2_output.fasta")
    copyfile(contig_out_file, phase1_2_output_file)
    
    # Process long read sequences
    
    lr_out_file = os.path.join(output_dir, "longread_IDS.fasta")
    with open(lr_out_file, 'w') as lr_fasta_file:
        with open(lr_fasta) as lf:
            write_flag = False
            for line in lf:
                line = line.strip()
                if line.startswith(">"):
                    header = line.split()[0][1:]
                    write_flag = header in lr_ids
                    if write_flag:
                        lr_fasta_file.write(line + "\n")
                elif write_flag:
                    lr_fasta_file.write(line + "\n")
    
    # Copy lr_out_file to unusedlongreads.fasta
    unused_longreads_file = os.path.join(input_dir, "unusedlongreads.fasta")
    copyfile(lr_out_file, unused_longreads_file)

def process_contig_ids(contig_set, input_dir, contig_fasta, lr_fasta, output_dir, folder_num):
    txt_files = set()
    for i in range(len(contig_set) - 1):
        txt_files.add(f"{contig_set[i]}_{contig_set[i+1]}.txt")
        txt_files.add(f"{contig_set[i+1]}_{contig_set[i]}.txt")
    
    lr_ids = set()
    for txt_file in txt_files:
        txt_file_path = os.path.join(input_dir, txt_file)
        if os.path.exists(txt_file_path):
            with open(txt_file_path) as f:
                for line in f:
                    c1, c2, lr_id = line.strip().split("\t")
                    if c1 in contig_set and c2 in contig_set:
                        lr_ids.add(lr_id)
    
    create_fasta_files(contig_set, lr_ids, contig_fasta, lr_fasta, output_dir, folder_num)

def main(input_file, input_dir, contig_fasta, lr_fasta, output_dir, num_processes, num_contigs):
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        folder_num = 1
        
        with open(input_file) as file:
            contig_ids = []
            for line in file:
                contig_ids.extend(line.strip().split(","))
                while len(contig_ids) >= num_contigs:
                    contig_set = contig_ids[:num_contigs]
                    futures.append(executor.submit(process_contig_ids, contig_set, input_dir, contig_fasta, lr_fasta, output_dir, folder_num))
                    contig_ids = contig_ids[num_contigs:]
                    folder_num += 1
            
            if contig_ids:
                futures.append(executor.submit(process_contig_ids, contig_ids, input_dir, contig_fasta, lr_fasta, output_dir, folder_num))
        
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error processing contig IDs: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 CreateBatchOfContigs.py input_file input_dir contig_fasta lr_fasta output_dir num_processes num_contigs")
        sys.exit(1)

    input_file = sys.argv[1]
    input_dir = sys.argv[2]
    contig_fasta = sys.argv[3]
    lr_fasta = sys.argv[4]
    output_dir = sys.argv[5]
    num_processes = int(sys.argv[6])
    num_contigs = int(sys.argv[7])

    main(input_file, input_dir, contig_fasta, lr_fasta, output_dir, num_processes, num_contigs)

print("Done!!!")

'''
import os
import sys
import concurrent.futures
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def create_fasta_files(contig_ids, lr_ids, contig_fasta, lr_fasta, output_dir, folder_num):
    folder_name = str(folder_num)
    output_dir = os.path.join(output_dir, folder_name)
    os.makedirs(output_dir, exist_ok=True)

    
    contig_out_file = os.path.join(output_dir, f"{folder_name}_Contigs.fasta")
    contig_records = SeqIO.to_dict(SeqIO.parse(contig_fasta, "fasta"))
    contig_seqs = []
    for contig_id in contig_ids:
        if contig_id in contig_records:
            contig_seq = contig_records[contig_id].seq
            contig_seqs.append(SeqRecord(contig_seq, id=contig_id, description=""))
        else:
            print(f"Contig {contig_id} not found in the contig FASTA file.")
            
    with open(contig_out_file, 'w') as contig_fasta_file:
        SeqIO.write(contig_seqs, contig_fasta_file, "fasta")

    lr_out_file = os.path.join(output_dir, "longread_IDS.fasta")
    with open(lr_out_file, 'w') as lr_fasta_file:
        with open(lr_fasta) as lf:
            write_flag = False
            for line in lf:
                line = line.strip()
                if line.startswith(">"):
                    header = line.split()[0][1:]
                    if header in lr_ids:
                        write_flag = True
                        lr_fasta_file.write(line + "\n")
                    else:
                        write_flag = False
                elif write_flag:
                    lr_fasta_file.write(line + "\n")

def process_contig_ids(contig_ids, num_contigs, input_dir, contig_fasta, lr_fasta, output_dir, folder_num):
    contig_sets = []
    for i in range(0, len(contig_ids), num_contigs):
        contig_set = contig_ids[i:i + num_contigs]
        contig_sets.append(contig_set)

    for i, contig_set in enumerate(contig_sets):
        if i < len(contig_sets) - 1:
            txt_files = set()
            for j in range(len(contig_set) - 1):
                contig_id1 = contig_set[j]
                contig_id2 = contig_set[j + 1]
                txt_file1 = f"{contig_id1}_{contig_id2}.txt"
                txt_file2 = f"{contig_id2}_{contig_id1}.txt"
                txt_files.add(txt_file1)
                txt_files.add(txt_file2)
                
            lr_ids = set()
            for txt_file in txt_files:
                txt_file_path = os.path.join(input_dir, txt_file)
                if not os.path.exists(txt_file_path):
                    continue
                with open(txt_file_path) as f:
                    for line in f:
                        c1, c2, lr_id = line.strip().split("\t")
                        if c1 in contig_set and c2 in contig_set:
                            lr_ids.add(lr_id)

            create_fasta_files(contig_set, lr_ids, contig_fasta, lr_fasta, output_dir, folder_num)
        else:
            remaining_contigs = num_contigs - len(contig_set)
            next_row_contig_ids = contig_ids[num_contigs:]
            if len(next_row_contig_ids) >= remaining_contigs:
                contig_set.extend(next_row_contig_ids[:remaining_contigs])
                contig_ids = next_row_contig_ids[remaining_contigs:]
            else:
                contig_set.extend(next_row_contig_ids)
                contig_ids = []

            txt_files = set()
            for j in range(len(contig_set) - 1):
                contig_id1 = contig_set[j]
                contig_id2 = contig_set[j + 1]
                txt_file1 = f"{contig_id1}_{contig_id2}.txt"
                txt_file2 = f"{contig_id2}_{contig_id1}.txt"
                txt_files.add(txt_file1)
                txt_files.add(txt_file2)

            lr_ids = set()
            for txt_file in txt_files:
                txt_file_path = os.path.join(input_dir, txt_file)
                if not os.path.exists(txt_file_path):
                    continue
                with open(txt_file_path) as f:
                    for line in f:
                        c1, c2, lr_id = line.strip().split("\t")
                        if c1 in contig_set and c2 in contig_set:
                            lr_ids.add(lr_id)

            create_fasta_files(contig_set, lr_ids, contig_fasta, lr_fasta, output_dir, folder_num)

def main(input_file, input_dir, contig_fasta, lr_fasta, output_dir, num_processes, num_contigs):
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = []

        with open(input_file) as input_file:
            contig_ids = []
            folder_num = 1
            for line in input_file:
                contig_ids.extend(line.strip().split(","))
                while len(contig_ids) >= num_contigs:
                    contig_set = contig_ids[:num_contigs]
                    futures.append(executor.submit(process_contig_ids, contig_set, num_contigs, input_dir, contig_fasta, lr_fasta, output_dir, folder_num))
                    contig_ids = contig_ids[num_contigs:]
                    folder_num += 1

        if contig_ids:
            futures.append(executor.submit(process_contig_ids, contig_ids, num_contigs, input_dir, contig_fasta, lr_fasta, output_dir, folder_num))

        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error processing contig IDs: {e}")

###LRs To Contigs Jem Wiring
if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python3 CreateBatchOfContigs.py input_file input_dir contig_fasta lr_fasta output_dir num_processes num_contigs")
        sys.exit(1)

    input_file = sys.argv[1]
    input_dir = sys.argv[2]
    contig_fasta = sys.argv[3]
    lr_fasta = sys.argv[4]
    output_dir = sys.argv[5]
    num_processes = int(sys.argv[6])
    num_contigs = int(sys.argv[7])

    main(input_file, input_dir, contig_fasta, lr_fasta, output_dir, num_processes, num_contigs)

print("Done!!!")
'''
