#!/bin/bash

# Help function
print_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -c, --contigs    Path to the contigs input file"
    echo "  -lr, --longreads Path to the long reads input file"
    echo "  -o, --output     Output directory (default: \$HOME/Maptcha/Output/)"
    echo "  -t, --threads    Number of threads to use (default: 16)"
    echo "  -n, --nodes      Number of nodes to use (default: 2)"
    echo "  -p, --processes  Number of processes per node (default: 2)"
    echo "  -h, --help       Show this help message"
}

# Default values
output_dir="$HOME/Maptcha/Output/"
threads=16
nodes=2
processes=2

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--contigs)
            contigs_input_file="$2"
            shift 2
            ;;
        -lr|--longreads)
            long_reads_input_file="$2"
            shift 2
            ;;
        -o|--output)
            output_dir="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -n|--nodes)
            nodes="$2"
            shift 2
            ;;
        -p|--processes)
            processes="$2"
            shift 2
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_help
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$contigs_input_file" || -z "$long_reads_input_file" ]]; then
    echo "Error: contigs and long reads input files are required."
    print_help
    exit 1
fi

# Check if the input files exist
if [[ ! -f "$contigs_input_file" ]]; then
    echo "Contigs input file '$contigs_input_file' not found."
    exit 1
fi

if [[ ! -f "$long_reads_input_file" ]]; then
    echo "Long reads input file '$long_reads_input_file' not found."
    exit 1
fi

# Recreate output directory
rm -rf "$output_dir"
mkdir -p "$output_dir"

# Detect cluster environment
if [[ ! -z "$SLURM_JOB_ID" ]]; then
    cluster="SLURM"
elif command -v qsub &> /dev/null; then
    cluster="PBS"
else
    echo "Unsupported cluster environment. Neither SLURM nor PBS commands found."
    exit 1
fi

# Calculate the total number of processes
np=$((nodes * processes))
# Main processing steps
chmod +x $HOME/Maptcha/src/CreateFastaFromLR
chmod +x $HOME/Maptcha/src/jem
$HOME/Maptcha/src/CreateFastaFromLR "$long_reads_input_file" "$output_dir/lr_leftright.fa" "$output_dir/lr_concat.fa"

start_time_maptcha=$(date +%s)
#chmod +x $HOME/Maptcha/src/CreateFastaFromLR
#chmod +x $HOME/Maptcha/src/jem
#$HOME/Maptcha/src/CreateFastaFromLR "$long_reads_input_file" "$output_dir/lr_leftright.fa" "$output_dir/lr_concat.fa"
mpiexec -np $np $HOME/Maptcha/src/jem -s "$contigs_input_file" -q "$output_dir/lr_concat.fa" -a $HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/A.txt -b /$HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/B.txt -p $HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/Prime.txt -r 1000 -n 30
cd ~/Maptcha/TestInput/
map_output="$HOME/Maptcha/TestInput/CLPairs.log"

$HOME/Maptcha/bin/GraphConstrWH "$map_output" "$output_dir/graphWH.txt"
$HOME/Maptcha/bin/graphLRID "$map_output" "$output_dir/graphLRID.txt"

python3 $HOME/Maptcha/src/Wiring.py "$map_output" "$output_dir/graphWH.txt" "$output_dir/wired_output.txt"
python3 $HOME/Maptcha/src/PathEnumeration.py "$map_output" "$output_dir/wired_output.txt" "$output_dir/path_output.txt"
python3 $HOME/Maptcha/src/CreateCCFileForEachPair_args.py "$output_dir/wired_output.txt" "$output_dir/graphLRID.txt" "$output_dir"
#python3 $HOME/Maptcha/src/CreateBatchOfContigs.py "$output_dir/path_output.txt" "$output_dir" "$contigs_input_file" "$long_reads_input_file" "$output_dir/FastaFilesBatch_8192/" 40 8192
python3 $HOME/Maptcha/src/CreateBatchOfContigs.py "$output_dir/path_output.txt" "$output_dir" "$contigs_input_file" "$long_reads_input_file" "$output_dir/FastaFilesBatch_8192/" $processes 8192

# Job submission
cd $HOME/Maptcha/Hifiasm/
chmod +x $HOME/Maptcha/Hifiasm/hifiasm

input_dir="$output_dir/FastaFilesBatch_8192"
mkdir -p "$output_dir/jobScripts/"
job_scripts_dir="$output_dir/jobScripts/"

if [[ -d "$job_scripts_dir" ]]; then
    echo "Job scripts directory created successfully."
else
    echo "Failed to create job scripts directory."
    exit 1
fi

start_time=$(date +%s)
all_folders=($(ls -d "$input_dir"/*/))
num_batches=${#all_folders[@]}

batches_per_node=$(( (num_batches + nodes - 1) / nodes ))
time_list=()

calculate_stats() {
    local sum_time=0
    local num_folders=${#time_list[@]}
    if ((num_folders > 0)); then
        for folder_time in "${time_list[@]}"; do
            sum_time=$((sum_time + folder_time))
        done
        average_time=$((sum_time / num_folders))

        if ((num_folders > 1)); then
            sum_squared_deviations=0
            for folder_time in "${time_list[@]}"; do
                deviation=$((folder_time - average_time))
                sum_squared_deviations=$((sum_squared_deviations + (deviation * deviation)))
            done
            variance=$((sum_squared_deviations / num_folders))
            standard_deviation=$(printf "%.2f" "$(echo "scale=10; sqrt($variance)" | bc)")
        fi
    fi
}

for ((node=0; node < nodes; node++)); do
    start=$((node * batches_per_node))
    end=$(( (node + 1) * batches_per_node ))
    if ((end >= num_batches)); then
        end=$((num_batches))
    fi

    batch_folders=("${all_folders[@]:start:end}")
    job_script="${job_scripts_dir}/job_script_node_${node}.sh"

    cat > "$job_script" <<EOF
#!/bin/bash
EOF

    if [[ "$cluster" == "PBS" ]]; then
        cat >> "$job_script" <<EOF
#PBS -l nodes=1:ppn=$processes
#PBS -l mem=120gb
#PBS -l walltime=06:00:00

cd $HOME/Maptcha/Hifiasm
chmod +x $HOME/Maptcha/Hifiasm/hifiasm

EOF
    elif [[ "$cluster" == "SLURM" ]]; then
        cat >> "$job_script" <<EOF
#SBATCH --nodes=1
#SBATCH --ntasks=$processes
#SBATCH --mem=120G
#SBATCH --time=06:00:00

cd $HOME/Maptcha/Hifiasm
chmod +x $HOME/Maptcha/Hifiasm/hifiasm

EOF
    fi

    for folder in "${batch_folders[@]}"; do
        folder_name=$(basename "$folder")
        base_folder_name="${folder_name%%.*}"
        contigs_file="$folder/${folder_name}_Contigs.fasta"
        long_reads_file="$folder/longread_IDS.fasta"
        output_file="${folder_name}.asm"

        cat >> "$job_script" <<EOF
start_folder_time=\$(date +%s)
#$HOME/Maptcha/Hifiasm/hifiasm -o "$folder/$output_file" -t $threads -n1 -a1 -r1 -f0 "$contigs_file" "$long_reads_file" > /dev/null 2>&1
end_folder_time=\$(date +%s)
folder_time=\$((end_folder_time - start_folder_time))
echo "Time taken for folder ${folder_name}: \${folder_time} seconds"
EOF
    done

    cat >> "$job_script" <<EOF
end_time=\$(date +%s)
total_time=\$((end_time - start_time))
echo "Total time taken: \${total_time} seconds"
EOF

    if [[ "$cluster" == "PBS" ]]; then
        qsub "$job_script"
    elif [[ "$cluster" == "SLURM" ]]; then
        sbatch "$job_script"
    fi
done

chmod +x $HOME/Maptcha/src/maptcha
$HOME/Maptcha/src/maptcha "$contigs_input_file" "$long_reads_input_file" "$output_dir"
calculate_stats
echo "Batched assembly done! "
# Calculate the total elapsed time for creating and submitting all job scripts
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

# Print the total time taken
echo "Job scripts created and submitted in $elapsed_time seconds."
$HOME/Maptcha/src/CreateFastaFromLR "$long_reads_input_file" "$output_dir/lr_leftright.fa" "$output_dir/lr_concat.fa"
mpiexec -np $np $HOME/Maptcha/src/jem -s "$output_dir/contExp.fasta" -q "$output_dir/lr_concat.fa" -a $HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/A.txt -b /$HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/B.txt -p $HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/Prime.txt -r 1000 -n 30
python3 $HOME/Maptcha/src/CreateUnmappedUnusedLR.py "$output_dir/FastaFilesBatch_8192/" "$output_dir/contExp.fasta" "$long_reads_input_file" "$output_dir/unused_longreads.fasta"
chmod +x $HOME/Maptcha/src/jem
cd $HOME/Maptcha/Hifiasm/
chmod +x $HOME/Maptcha/Hifiasm/hifiasm

# Start the timer
start_time=$(date +%s)
mkdir -p "$output_dir/Phase2/"
$HOME/Maptcha/Hifiasm/hifiasm -o "$output_dir/Phase2/Only_UnmappedUnusedLongreads.asm" -t $threads -n1 -a1 -r1 -f0 "$output_dir/unused_longreads.fasta" > /dev/null 2>&1
awk '/^S/{print ">"$2;print $3}' "$output_dir/Phase2/Only_UnmappedUnusedLongreads.asm.bp.p_ctg.gfa" > "$output_dir/Phase2/Only_UnmappedUnusedLongreads.asm.bp.p_ctg.gfa.fa"
# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Longread Island Construction done! "
#rm "$output_dir/contExp.fasta"
python3 $HOME/Maptcha/src/merge.py "$output_dir/Phase2/Only_UnmappedUnusedLongreads.asm.bp.p_ctg.gfa.fa" "$output_dir/contExp.fasta" "$output_dir/Phase1_2_partialScaff.fa"
$HOME/Maptcha/src/CreateFastaFromLR "$l$output_dir/unused_longreads.fasta" "$output_dir/lr_leftright.fa" "$output_dir/lr_concat.fa"
rm "$output_dir/unusedlongreads.fasta"
rm "$output_dir/contExp.fasta"
mpiexec -np $np $HOME/Maptcha/src/jem -s "$output_dir/Phase1_2_partialScaff.fa" -q "$output_dir/lr_concat.fa" -a $HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/A.txt -b /$HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/B.txt -p $HOME/Maptcha/JEM-Mapper/TestInput/ConstantsForLCH/Prime.txt -r 1000 -n 30
rm "$output_dir/lr_leftright.fa" 
rm "$output_dir/lr_concat.fa"
cd $HOME/Maptcha/Hifiasm/
chmod +x $HOME/Maptcha/Hifiasm/hifiasm
# Start the timer
start_time=$(date +%s)
mkdir -p "$output_dir/Final/"
$HOME/Maptcha/Hifiasm/hifiasm -o "$output_dir/Final/finalAssembly.asm" -t $threads -n1 -a1 -r1 -f0 "$output_dir/unused_longreads.fasta" > /dev/null 2>&1
awk '/^S/{print ">"$2;print $3}' "$output_dir/Final/finalAssembly.asm.bp.p_ctg.gfa" > "$output_dir/Final/finalAssembly.fa"
rm "$output_dir/phase1_2_output.fasta"
rm "$output_dir/unused_longreads.fasta"
rm "$output_dir/Phase1_2_partialScaff.fa"
rm -rf "$output_dir/FastaFilesBatch_8192/"
rm -rf "$output_dir/jobScripts/"
calculate_stats() {
    local sum_time=0
    local num_folders=${#time_list[@]}
    if ((num_folders > 0)); then
        for folder_time in "${time_list[@]}"; do
            sum_time=$((sum_time + folder_time))
        done
        average_time=$((sum_time / num_folders))

        if ((num_folders > 1)); then
            sum_squared_deviations=0
            for folder_time in "${time_list[@]}"; do
                deviation=$((folder_time - average_time))
                sum_squared_deviations=$((sum_squared_deviations + (deviation * deviation)))
            done
            variance=$((sum_squared_deviations / num_folders))
            standard_deviation=$(printf "%.2f" "$(echo "scale=10; sqrt($variance)" | bc)")
        fi
    fi
}
# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Link Scaffolds With Bridges done! "

end_time_maptcha=$(date +%s)
elapsed_time_maptcha=$((end_time_maptcha - start_time_maptcha))
#echo "The total time is: $elapsed_time_maptcha seconds using $threads thread, $nodes nodes and $processes processes."
echo "The total time is: $elapsed_time_maptcha seconds "
echo "Thank you for using Maptcha! "
