// Maptcha: An efficient parallel workflow for hybrid genome scaffolding

// Oieswarya Bhowmik, Tazin Rahman, Ananth Kalyanaraman

//      (oieswarya.bhowmik@wsu.edu, tazin.rahman@wsu.edu, ananth@wsu.edu)

// Washington State University

//

// **************************************************************************************************

// Copyright (c) 2024. Washington State University ("WSU"). All Rights Reserved.
// Permission to use, copy, modify, and distribute this software and its documentation
// for educational, research, and not-for-profit purposes, without fee, is hereby
// granted, provided that the above copyright notice, this paragraph and the following
// two paragraphs appear in all copies, modifications, and distributions. For
// commercial licensing opportunities, please contact The Office of Commercialization,
// WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
// commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

// IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
// THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
// ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
// OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
// **************************************************************************************************



import os
import sys

def create_wiring_files(file1, file2, output_dir):
    # read in File1 and create a dictionary mapping contig pairs to LRs
    contig_pairs = {}

    with open(file1, "r") as f1:
        next(f1)  # skip first line
        next(f1)  # skip second line
        for line in f1:
            cols = line.strip().split()
            contig1, contig2, lr_count = cols[0], cols[1], cols[2]
            contig_pair = tuple(sorted([contig1, contig2]))
            contig_pairs[contig_pair] = int(lr_count)

    # read in File2 and create a dictionary mapping contig pairs to LR IDs
    lr_ids = {}

    with open(file2, "r") as f2:
        for line in f2:
            cols = line.strip().split()
            contig1, contig2, lr_id = cols[0], cols[1], cols[2]
            contig_pair = tuple(sorted([contig1, contig2]))
            if contig_pair in contig_pairs:
                if contig_pair not in lr_ids:
                    lr_ids[contig_pair] = set()
                lr_ids[contig_pair].add(lr_id)

    # create output directory if it doesn't already exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # write output files for each contig pair
    for contig_pair, lr_id_set in lr_ids.items():
        contig1, contig2 = contig_pair
        out_file = os.path.join(output_dir, f"{contig1}_{contig2}.txt")
        with open(out_file, "w") as f_out:
            for lr_id in lr_id_set:
                f_out.write(f"{contig1}\t{contig2}\t{lr_id}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 /home/obhowmik/PythonPractice/CreateCCFileForEachPair_args.py file1 file2 output_directory")
        sys.exit(1)

    create_wiring_files(sys.argv[1], sys.argv[2], sys.argv[3])
