#!/usr/bin/env python3
import sys
import os
import json
import argparse
from pyfaidx import Fasta
from multiprocessing import Pool
import py_kmc_api as pka
import traceback

# Global constants 
FLANKING_NUCLEOTIDES = 30
NUM_PROCESSES = 20
PROGRESS_INTERVAL = 5000

# Global variables that will be set via CLI arguments
K = None
DATA_DIR = None
MAF_FILEPATH = None
FASTA_FILE = None
NEOMERS_OUTPUT_DIR = None
kmer_data_base = None


def checkIfExistsInKmerSet(kmer_seq):
    global kmer_data_base
    kmer_obj = pka.KmerAPI(K)
    kmer_obj.from_string(kmer_seq)
    exists = kmer_data_base.IsKmer(kmer_obj)
    return exists

def apply_mutation(original_ref_seq, variant_type, ref_allele, alt_allele, debug_msgs):
    print("\n[apply_mutation] --- START ---")
    print(f"[apply_mutation] Variant Type: {variant_type}")
    print(f"[apply_mutation] Reference Allele (from MAF): '{ref_allele}'")
    print(f"[apply_mutation] Alternate Allele (from MAF / Tumor_Seq_Allele2): '{alt_allele}'")
    print(f"[apply_mutation] Original Reference Sequence (from genome): '{original_ref_seq}'")
    print("[apply_mutation] The original_ref_seq should start with the ref_allele if ref_allele != '-'.")

    # Validate ref_allele
    if ref_allele != "-" and not original_ref_seq.startswith(ref_allele):
        error_str = (f"[apply_mutation] ERROR: Variant Type: {variant_type} Applying Tumor_Seq_Allele2: '{alt_allele}'  "
                     f"Reference allele mismatch! Expected '{ref_allele}' at start of '{original_ref_seq}'")
        print(error_str)
        raise ValueError(f"Reference allele does not match the original reference sequence from genome.\n{error_str}\n{debug_msgs}")

    # Apply logic based on variant type
    if variant_type in ("SNP", "DNP", "TNP", "ONP"):
        if len(ref_allele) != len(alt_allele):
            print(f"[apply_mutation] ERROR: Length mismatch in substitution: ref='{ref_allele}', alt='{alt_allele}'")
            raise ValueError("Substitution allele length mismatch.")
        mutated_seq = alt_allele + original_ref_seq[len(ref_allele):]
        print("[apply_mutation] Substitution applied. "
              f"Replacing '{ref_allele}' with '{alt_allele}' results in '{mutated_seq}'.")

    elif variant_type == "INS":
        if ref_allele != "-":
            print(f"[apply_mutation] ERROR: For insertion expected ref_allele='-', got '{ref_allele}'.")
            raise ValueError("Insertion requires ref_allele='-'.")
        mutated_seq = alt_allele + original_ref_seq
        print("[apply_mutation] Insertion applied. "
              f"Inserting '{alt_allele}' before '{original_ref_seq}' results in '{mutated_seq}'.")

    elif variant_type == "DEL":
        if alt_allele != "-":
            print(f"[apply_mutation] ERROR: For deletion expected alt_allele='-', got '{alt_allele}'.")
            raise ValueError("Deletion requires alt_allele='-'.")
        mutated_seq = original_ref_seq[len(ref_allele):]
        print("[apply_mutation] Deletion applied. "
              f"Removing '{ref_allele}' from '{original_ref_seq}' results in '{mutated_seq}'.")

    else:
        print(f"[apply_mutation] WARNING: Unsupported variant type '{variant_type}'. Returning original sequence.")
        mutated_seq = original_ref_seq

    print(f"[apply_mutation] Resulting Mutated Sequence: '{mutated_seq}'")
    print("[apply_mutation] --- END ---\n")
    return mutated_seq

def get_kmers(seq, k):
    length = len(seq)
    print(f"[get_kmers] Generating {k}-mers from sequence of length {length}.")
    kmers = {seq[i:i+k] for i in range(length - k + 1)} if length >= k else set()
    print(f"[get_kmers] Generated {len(kmers)} k-mers.")
    return kmers

def reverse_complement(seq):
    """
    Generates the reverse complement of a DNA sequence.
    """
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def process_worker(args_tuple):
    # Each worker now initializes its own fasta handle:
    local_fasta = Fasta(f'{DATA_DIR}/{FASTA_FILE}')
    global K, FLANKING_NUCLEOTIDES

    (lines, header_indices, worker_id, total_lines, output_dir) = args_tuple
    (idx_Hugo_Symbol, idx_Chromosome, idx_Start_position, idx_End_position, idx_Strand,
     idx_Variant_Classification, idx_Variant_Type, idx_Reference_Allele, idx_Tumor_Seq_Allele1,
     idx_Tumor_Seq_Allele2, idx_Tumor_Sample_Barcode, idx_Matched_Norm_Sample_Barcode,
     idx_Donor_ID, idx_Project_Code) = header_indices

    output_file = os.path.join(output_dir, f"worker_{worker_id}.txt")
    error_file = os.path.join(output_dir, f"error_worker_{worker_id}.err")
    print(f"[process_worker] Worker {worker_id}: Start processing. Output file: {output_file}, Error file: {error_file}")

    with open(output_file, "w") as out, open(error_file, "w") as err_out:
        processed = 0
        potential_neomers_count = 0
        final_neomers_count = 0
        final_neomers_set = set()
        total_to_process = len(lines)

        for line in lines:
            cols = line.strip().split('\t')
            try:
                Hugo_Symbol = cols[idx_Hugo_Symbol]
                Chromosome = cols[idx_Chromosome]
                Start_position = int(cols[idx_Start_position]) - 1
                End_position = int(cols[idx_End_position])
                Strand = cols[idx_Strand]
                Variant_Classification = cols[idx_Variant_Classification]
                Variant_Type = cols[idx_Variant_Type]
                Reference_Allele = cols[idx_Reference_Allele]
                Tumor_Seq_Allele1 = cols[idx_Tumor_Seq_Allele1]
                Tumor_Seq_Allele2 = cols[idx_Tumor_Seq_Allele2]
                Tumor_Sample_Barcode = cols[idx_Tumor_Sample_Barcode]
                Matched_Norm_Sample_Barcode = cols[idx_Matched_Norm_Sample_Barcode]
                Donor_ID = cols[idx_Donor_ID]
                Project_Code = cols[idx_Project_Code]

                if not Chromosome.startswith("chr"):
                    Chromosome = "chr" + Chromosome

                print(f"\n[process_worker] Worker {worker_id}: Processing mutation:")
                print(f"  Chromosome: {Chromosome}")
                print(f"  Start_position (0-based): {Start_position}")
                print(f"  End_position (1-based): {End_position}")
                print(f"  Variant_Type: {Variant_Type}")
                print(f"  Reference_Allele: '{Reference_Allele}', Tumor_Seq_Allele1: '{Tumor_Seq_Allele1}', Tumor_Seq_Allele2: '{Tumor_Seq_Allele2}'")

                alt_allele = Tumor_Seq_Allele2
                print(f"[process_worker] Determined ALT allele is always Tumor_Seq_Allele2: '{alt_allele}'")

                # Extract original reference sequence
                original_ref_seq = local_fasta[Chromosome][Start_position:End_position].seq.upper()
                if len(original_ref_seq) == 0:
                    error_msg = f"[ERROR] Empty original_ref_seq for {Chromosome}:{Start_position}-{End_position}. Check coordinates!\n"
                    print(error_msg)
                    raise ValueError("No reference bases extracted from genome. Coordinates may be invalid.")

                orig_flank_left = max(Start_position - 5, 0)
                orig_flank_right = min(End_position + 5, len(local_fasta[Chromosome]))
                orig_flank_seq = local_fasta[Chromosome][orig_flank_left:orig_flank_right].seq.upper()
                debug_msg1 = f"[process_worker] Original reference flank ({Chromosome}:{orig_flank_left}-{orig_flank_right}): '{orig_flank_seq}'"
                debug_msg2 = f"[process_worker] Extracted original_ref_seq: '{original_ref_seq}' (Expected to start with '{Reference_Allele}' if Ref_Allele != '-')"
                print(debug_msg1)
                print(debug_msg2)

                mutated_seq = apply_mutation(original_ref_seq, Variant_Type, Reference_Allele, alt_allele,
                                             f"debug1: {debug_msg1}\ndebug2: {debug_msg2}\n")

                chrom_len = len(local_fasta[Chromosome])
                FLANKING = K + FLANKING_NUCLEOTIDES  # = (K + 30)
                window_start = max(0, Start_position - FLANKING)
                window_end = min(chrom_len, End_position + FLANKING)
                window_seq = local_fasta[Chromosome][window_start:window_end].seq.upper()

                mut_start_idx = Start_position - window_start
                original_mutation_length = len(original_ref_seq)

                mutated_window_seq = (window_seq[:mut_start_idx] +
                                      mutated_seq +
                                      window_seq[mut_start_idx + original_mutation_length:])

                print(f"[process_worker] Window coordinates: {Chromosome}:{window_start}-{window_end}")
                print(f"[process_worker] Window seq length: {len(window_seq)}")
                print(f"[process_worker] Window_Seq='{window_seq}'")
                print(f"[process_worker] Mutated_Window_Seq='{mutated_window_seq}'")
                print(f"[process_worker] Mutation start in window: {mut_start_idx}, end: {mut_start_idx + len(mutated_seq)}")

                check_flank_size = 5
                check_start = max(0, mut_start_idx - check_flank_size)
                check_end = min(len(window_seq), mut_start_idx + len(mutated_seq) + check_flank_size)
                original_check_seq = window_seq[check_start:check_end]
                mutated_check_seq = mutated_window_seq[check_start:check_end]
                print(f"[process_worker] Original check seq around mutation: '{original_check_seq}'")
                print(f"[process_worker] Mutated check seq around mutation:  '{mutated_check_seq}'")

                original_kmers = get_kmers(window_seq, K)
                mutated_kmers = get_kmers(mutated_window_seq, K)

                potential_neomers = mutated_kmers - original_kmers
                potential_neomers_count += len(potential_neomers)
                print(f"[process_worker] Potential neomers: {len(potential_neomers)}")
                
                # Append Reverse Complements to potential_neomers
                rev_complements = [reverse_complement(kmer) for kmer in potential_neomers]
                potential_neomers = list(potential_neomers)
                potential_neomers += rev_complements
                print(f"[process_worker] Worker {worker_id}: Appended reverse complements to potential_neomers.")
                print(f"[process_worker] Potential neomers count before filtering: {len(potential_neomers) - len(rev_complements)}")
                print(f"[process_worker] Potential neomers count after appending reverse complements: {len(potential_neomers)}")
                
                # Now filter the augmented potential_neomers
                actual_neomers = [nm for nm in potential_neomers if not checkIfExistsInKmerSet(nm)]
                final_neomers_count += len(actual_neomers)
                final_neomers_set.update(actual_neomers)
                print(f"[process_worker] Worker {worker_id}: Identified {len(actual_neomers)} actual neomers after filtering.")

                if actual_neomers:
                    neomers_str = ",".join(actual_neomers)
                    out.write(f"{Donor_ID}\t{Project_Code}\t{Hugo_Symbol}\t{Variant_Classification}\t{neomers_str}\n")

                # Every progress interval print completion percentage
                processed += 1
                if processed % PROGRESS_INTERVAL == 0:
                    percent = (processed / total_to_process) * 100
                    print(f"[process_worker] Worker {worker_id}: {processed}/{total_to_process} processed "
                          f"({percent:.2f}%). Potential Neomers: {potential_neomers_count}, Final Neomers: {final_neomers_count}")

            except Exception as e:
                # This 'except' corresponds to the 'try' inside the for loop
                err_out.write(f"Error processing line: {line}\n")
                err_out.write(f"Exception: {str(e)}\n")
                err_out.write(traceback.format_exc())
                err_out.write("\n")
                continue  # Continue with the next line

    print(f"[process_worker] Worker {worker_id} completed.")
    print(f"[process_worker] Processed: {processed}")
    print(f"[process_worker] Potential Neomers: {potential_neomers_count}")
    print(f"[process_worker] Final Neomers: {final_neomers_count}")
    print(f"[process_worker] Unique Final Neomers Found: {len(final_neomers_set)}")
    return (worker_id, output_file, total_to_process, processed, potential_neomers_count, final_neomers_count, final_neomers_set)

def main(args):
    global K, DATA_DIR, MAF_FILEPATH, FASTA_FILE, NEOMERS_OUTPUT_DIR, kmer_data_base

    # Assign CLI argument values to global variables
    K = args.K
    DATA_DIR = args.data_dir
    MAF_FILEPATH = args.maf_filepath
    FASTA_FILE = args.fasta_file
    NEOMERS_OUTPUT_DIR = args.neomers_output_dir
    OUTFILE = args.outfile
    if OUTFILE == './K_neomers.output.txt':
        OUTFILE = f"{K}_neomers_output.txt"
    print("[main] Starting kmer processing...")
    print(f"[main] Parsing MAF file {MAF_FILEPATH} to detect header and indexing...")

    # Open the KMC database using the provided DATA_DIR and K
    KMC_DATABASE_PATH = f'{DATA_DIR}/{K}mers.res'
    kmer_data_base = pka.KMCFile()
    if not kmer_data_base.OpenForRA(KMC_DATABASE_PATH):
        print("[ERROR] Cannot open KMC database")
        sys.exit(1)

    header_line = None
    with open(MAF_FILEPATH, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('Hugo_Symbol'):
                header_line = line.strip()
                break

    if header_line is None:
        print("[main] ERROR: No header found in MAF file. Cannot proceed.")
        sys.exit(1)

    headers = header_line.split('\t')
    print(f"[main] Detected headers: {headers}")

    def col_idx(colname):
        try:
            return headers.index(colname)
        except ValueError:
            print(f"[main] ERROR: Column {colname} not found in MAF headers.")
            sys.exit(1)

    idx_Hugo_Symbol = col_idx('Hugo_Symbol')
    idx_Chromosome = col_idx('Chromosome')
    idx_Start_position = col_idx('Start_position')
    idx_End_position = col_idx('End_position')
    idx_Strand = col_idx('Strand')
    idx_Variant_Classification = col_idx('Variant_Classification')
    idx_Variant_Type = col_idx('Variant_Type')
    idx_Reference_Allele = col_idx('Reference_Allele')
    idx_Tumor_Seq_Allele1 = col_idx('Tumor_Seq_Allele1')
    idx_Tumor_Seq_Allele2 = col_idx('Tumor_Seq_Allele2')
    idx_Tumor_Sample_Barcode = col_idx('Tumor_Sample_Barcode')
    idx_Matched_Norm_Sample_Barcode = col_idx('Matched_Norm_Sample_Barcode')
    idx_Donor_ID = col_idx('Donor_ID')
    idx_Project_Code = col_idx('Project_Code')

    header_indices = (idx_Hugo_Symbol, idx_Chromosome, idx_Start_position, idx_End_position, idx_Strand,
                      idx_Variant_Classification, idx_Variant_Type, idx_Reference_Allele, idx_Tumor_Seq_Allele1,
                      idx_Tumor_Seq_Allele2, idx_Tumor_Sample_Barcode, idx_Matched_Norm_Sample_Barcode,
                      idx_Donor_ID, idx_Project_Code)

    all_lines = []
    with open(MAF_FILEPATH, 'r') as f:
        start_collecting = False
        for line in f:
            if start_collecting:
                if line.startswith('#') or line.startswith('Hugo_Symbol'):
                    continue
                all_lines.append(line)
            elif line.startswith('Hugo_Symbol'):
                start_collecting = True

    total_lines = len(all_lines)
    print(f"[main] Total MAF mutations (lines) to process: {total_lines}")

    output_dir = "worker_outputs"
    os.makedirs(output_dir, exist_ok=True)
    chunk_size = total_lines // NUM_PROCESSES
    chunks = []
    start = 0
    for wid in range(NUM_PROCESSES):
        end = start + chunk_size
        if wid == NUM_PROCESSES - 1:
            end = total_lines
        chunk_lines = all_lines[start:end]
        start = end
        chunks.append((chunk_lines, header_indices, wid, total_lines, output_dir))
        print(f"[main] Assigned {len(chunk_lines)} lines to Worker {wid}")

    print(f"[main] Using {NUM_PROCESSES} processes for multiprocessing.")
    with Pool(NUM_PROCESSES) as pool:
        results = pool.map(process_worker, chunks)

    if not os.path.exists(NEOMERS_OUTPUT_DIR):
        os.makedirs(NEOMERS_OUTPUT_DIR)
    output_file = f"{NEOMERS_OUTPUT_DIR}/{K}_neomers_output.txt"
    print(f"[main] Merging all worker output files into {output_file}")

    with open(output_file, 'w') as out:
        out.write("Donor_ID\tProject_Code\tHugo_Symbol\tVariant_Classification\tneomers_created\n")
        for (wid, worker_file, tot, proc, pnc, fnc, fns) in results:
            if os.path.getsize(worker_file) > 0:
                with open(worker_file, 'r') as tmp_in:
                    for line in tmp_in:
                        out.write(line)
            os.remove(worker_file)
            print(f"[main] Merged output from Worker {wid}: processed {proc} lines, Potential Neomers={pnc}, Final Neomers={fnc}")

    final_error_file = "error_mutations.err"
    with open(final_error_file, 'w') as final_err_out:
        for wid in range(NUM_PROCESSES):
            error_file = os.path.join(output_dir, f"error_worker_{wid}.err")
            if os.path.exists(error_file):
                with open(error_file, 'r') as ef:
                    for line in ef:
                        final_err_out.write(line)
                os.remove(error_file)
                print(f"[main] Merged error log from Worker {wid} into {final_error_file}")

    total_potential_neomers = sum(r[4] for r in results)
    total_final_neomers = sum(r[5] for r in results)
    all_neomers_sets = [r[6] for r in results]
    all_final_neomers = set()
    for s in all_neomers_sets:
        all_final_neomers.update(s)

    ratio = (total_final_neomers / total_potential_neomers * 100) if total_potential_neomers > 0 else 0

    print("== Final Statistics ==")
    print(f"Total Mutations: {total_lines}")
    print(f"Total Potential Neomers: {total_potential_neomers}")
    print(f"Total Actual Final Neomers: {total_final_neomers}")
    print(f"Percentage of potential neomers that are actual: {ratio:.2f}%")
    print(f"Total Unique Final Neomers Found: {len(all_final_neomers)}")

    for (wid, worker_file, tot, proc, pnc, fnc, fns) in results:
        percent = (proc / (tot / NUM_PROCESSES)) * 100 if tot > 0 else 100
        print(f"Worker {wid}: processed {proc} mutations (~{percent:.2f}% of assigned), "
              f"Potential neomers: {pnc}, Final neomers: {fnc}, Unique final neomers: {len(fns)}")

    print("[main] Processing completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Neomer CLI tool for processing mutations.")
    parser.add_argument('--K', type=int, default=11,
                        help='K value for k-mer generation (default: 11)')
    parser.add_argument('--data_dir', type=str, default='./data_working_dir',
                        help='Directory where data files (including KMC database and FASTA) are located (default: ./data_working_dir)')
    parser.add_argument('--maf_filepath', type=str,
                        default='./test_files/sample_1000_first_final_consensus_passonly.snv_mnv_indel.icgc.public.maf',
                        help='Path to the MAF file (default: ./test_files/sample_1000_first_final_consensus_passonly.snv_mnv_indel.icgc.public.maf)')
    parser.add_argument('--fasta_file', type=str, default='hg19.fa',
                        help='FASTA file name (default: hg19.fa)')
    parser.add_argument('--neomers_output_dir', type=str, default='./neomer_results',
                        help='Directory to write neomers output (default: ./neomer_results)')
    parser.add_argument('--outfile', type=str, default='./K_neomers.output.txt',
                        help='Default output filename (default: ./K_neomers.output.txt)')
    args = parser.parse_args()
    main(args)
