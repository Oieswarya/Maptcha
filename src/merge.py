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
