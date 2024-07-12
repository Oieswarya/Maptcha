# Maptcha: An efficient parallel workflow for hybrid genome scaffolding

# Oieswarya Bhowmik, Tazin Rahman, Ananth Kalyanaraman

#      (oieswarya.bhowmik@wsu.edu, tazin.rahman@wsu.edu, ananth@wsu.edu)

# Washington State University

#

# **************************************************************************************************

# Copyright (c) 2024. Washington State University ("WSU"). All Rights Reserved.
# Permission to use, copy, modify, and distribute this software and its documentation
# for educational, research, and not-for-profit purposes, without fee, is hereby
# granted, provided that the above copyright notice, this paragraph and the following
# two paragraphs appear in all copies, modifications, and distributions. For
# commercial licensing opportunities, please contact The Office of Commercialization,
# WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
# commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

# IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
# THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
# ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
# OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
# **************************************************************************************************



import sys

def merge_files(input_file1, input_file2, output_file):
    with open(output_file, "w") as output_f:
        with open(input_file1, "r") as file1:
            for line in file1:
                if line.startswith(">"):
                    header = line.strip()[1:]
                    modified_header = header + "_UnusedLRs"
                    output_f.write(">" + modified_header + "\n")
                else:
                    output_f.write(line)

        with open(input_file2, "r") as file2:
            for line in file2:
                if line.startswith(">"):
                    header = line.strip()[1:]
                    modified_header = header + "_conts"
                    output_f.write(">" + modified_header + "\n")
                else:
                    output_f.write(line)

def renumber_reads(input_file, output_file):
    with open(input_file, 'r') as data_file1:
        count = 0
        with open(output_file, "w") as output_file:
            for line in data_file1:
                line = line.strip()
                if line.startswith('>'):
                    count += 1
                    long_read_name = f">{count}"
                else:
                    long_read_name = line
                output_file.write(long_read_name + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merged_code.py input_file1 input_file2 output_file")
        sys.exit(1)

    input_file1 = sys.argv[1]
    input_file2 = sys.argv[2]
    output_file1 = sys.argv[3]

    intermediate_file = "intermediate_file.fa"

    merge_files(input_file1, input_file2, intermediate_file)
    renumber_reads(intermediate_file, output_file1)
    #print("Done")
