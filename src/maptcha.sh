#!/bin/bash
#PBS -l nodes=1:intel:ppn=1
#PBS -l mem=120gb
#PBS -l walltime=06:00:00

#Load the required modules
#module load gcc/7.3.0
#module load python/3.8/gcc/7.3.0
#module load sqlite/3.35.1/system
#module load openmpi/3.1.2/gcc/7.3.0


############      Creating the log file from JEM-Mapper output      ###############
python3 $HOME/Maptcha/src/CreateCLFromLog_PyCharm.py
cd ~/Maptcha/TestInput/
map_output="$HOME/Maptcha/TestInput/CLPairs.log"

#################################################################  Contig Expansion    ##############################################################################

################# Graph Construction  ###################

g++ -fopenmp -O3 -o ~/Maptcha/src/GraphConstrWH ~/Maptcha/src/GraphConstrWH.cpp
~/Maptcha/src/GraphConstrWH "$map_output" ~/Maptcha/graphWH.txt


g++ -fopenmp -O3 -o ~/Maptcha/src/graphLRID ~/Maptcha/src/graphLRID.cpp
~/Maptcha/src/graphLRID "$map_output" ~/Maptcha/graphLRID.txt


################# Wiring and Path Enumeration  ###################
 
python3 ~/Maptcha/src/Wiring.py "$map_output" ~/Maptcha/graphWH.txt "$HOME/Maptcha/wired_output.txt"


python3 ~/Maptcha/src/PathEnumeration.py "$map_output" $HOME/Maptcha/wired_output.txt $HOME/Maptcha/path_output.txt


################# Batching of Contigs  ###################

python3 $HOME/Maptcha/src/CreateCCFileForEachPair_args.py "$HOME/Maptcha/wired_output.txt" "$HOME/Maptcha/graphLRID.txt" "$HOME/Maptcha/Output/"


python3 $HOME/Maptcha/src/CreateBatchOfContigs.py "$HOME/Maptcha/path_output.txt" "$HOME/Maptcha/Output/" "$HOME/Maptcha/TestInput/minia_Coxiellaburnetii_contigs.fa" "$HOME/Maptcha//TestInput/CoxiellaBurnetii_longreads.fa" "$HOME/Maptcha/Output/FastaFilesBatch_8192/" 40 8192


## Define compute options for the main job
#PBS -l nodes=4:amd:ppn=2
#PBS -l mem=4gb
#PBS -l walltime=01:00:00

cd $HOME/Maptcha/Hifiasm/
chmod +x $HOME/Maptcha/Hifiasm/hifiasm



# Set the path to the directory containing the input folders
input_dir='$HOME/Maptcha/Output/FastaFilesBatch_8192/'

# Set the path to the job scripts directory
mkdir -p "$HOME/Maptcha/Output/jobScripts_Coxiella/"
job_scripts_dir="$HOME/Maptcha/Output/jobScripts_Coxiella/"

# Assign the job scripts directory to the variable
job_scripts_dir="$HOME/Maptcha/Output/jobScripts_Coxiella/"

# Check if the directory was created successfully
if [ -d "$job_scripts_dir" ]; then
    echo "Job scripts directory created successfully."
else
    echo "Failed to create job scripts directory."
    exit 1
fi

# Start the timer
start_time=$(date +%s)

# Get the list of all folders in the input directory
all_folders=($(ls -d "$input_dir"/*/))

# Count the total number of batches
num_batches=${#all_folders[@]}

# Get MPI rank
#MPI_RANK=$OMPI_COMM_WORLD_RANK
#MPI_SIZE=$OMPI_COMM_WORLD_SIZE

MPI_RANK=$(mpiexec -n 1 bash -c 'echo $OMPI_COMM_WORLD_RANK')

# Calculate the start and end indices of batches for this node
start_index=$((MPI_RANK * batches_per_node))
#end_index=$(( (MPI_RANK + 1) * batches_per_node ))
end_index=$(( start_index + batches_per_node - 1))


# Extract the batch folders for this node
batch_folders=("${all_folders[@]:start_index:end_index}")
  

# Set the number of nodes and processors per node
nodes=4
ppn=2

# Calculate the total number of processors available
np=$((nodes * ppn))

# Calculate the number of batches each node will handle
batches_per_node=$(( (num_batches + nodes - 1) / nodes ))

# Initialize variables for time calculations
total_time=0
max_time=0
time_list=()

# Function to calculate average and standard deviation
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

# Loop through each node
for ((node=0; node < nodes; node++)); do
  # Calculate the start and end indices of batches for this node
  start=$((node * batches_per_node))
  end=$(( (node + 1) * batches_per_node ))
  if ((end >= num_batches)); then
    end=$((num_batches))
  fi

  # Extract the batch folders for this node
  batch_folders=("${all_folders[@]:start:end}")

  # Create a unique job script for this node
  job_script="${job_scripts_dir}/job_script_node_${node}.sh"

  # Write the PBS directives and the script content to the job script
  cat > "$job_script" <<EOF
#!/bin/bash
## Define compute options for the node job
#PBS -l nodes=1:amd:ppn=$ppn
#PBS -l mem=120gb
#PBS -l walltime=06:00:00


cd $HOME/Maptcha/Hifiasm
chmod +x $HOME/Maptcha/Hifiasm/hifiasm


EOF

  # Loop through each batch and append the hifiasm command to the job script
  for folder in "${batch_folders[@]}"; do
    folder_name=$(basename "$folder")
    base_folder_name="${folder_name%%.*}" # Extract the base folder name without any file extension
    #contigs_file="$folder/${folder_name}_Moth.fasta"
    contigs_file="$folder/${folder_name}_Contigs.fasta"
    long_reads_file="$folder/longread_IDS.fasta"
    output_file="${folder_name}.asm" # Use the base folder name for the output file

    cat >> "$job_script" <<EOF
# Processing folder: $folder_name
start_folder_time=\$(date +%s)
$HOME/Maptcha/Hifiasm/hifiasm -o "$folder/$output_file" -t 1 "$contigs_file" "$long_reads_file"
end_folder_time=\$(date +%s)
folder_elapsed_time=\$((end_folder_time - start_folder_time))

# Output the time measurement to the error file for debugging
echo "Folder: $folder_name, Elapsed time: \$folder_elapsed_time seconds" >> "${job_script}.o"

awk '/^S/{print ">"\$2;print \$3}' "$folder/$output_file.bp.p_ctg.gfa" > "$folder/$output_file.bp.p_ctg.gfa.fa"
EOF
  done

  # Make the job script executable
  chmod +x "$job_script"

  # Submit the job script to the PBS scheduler
  qsub "$job_script" > "${job_script}.o"
done

# Wait for all batch jobs to finish
wait

# Initialize the time_list array
time_list=()
phase1_2_output="$HOME/Maptcha/TestInput/minia_Coxiellaburnetii_contigs.fa" 
unusedlongreads="$HOME/Maptcha/TestInput/CoxiellaBurnetii_longreads.fa"
# Loop through each batch and read the elapsed times from the output files
for ((node=0; node < nodes; node++)); do
  job_script="${job_scripts_dir}/job_script_node_${node}.sh"
  for ((i=start; i < end; i++)); do
    folder_name=$(basename "${all_folders[$i]}")
    elapsed_time=$(grep "Folder: $folder_name" "${job_script}.o" | awk '{print $NF}')
    time_list+=("$elapsed_time")
  done
done

# Calculate average and standard deviation
calculate_stats



echo "Batched assembly done! "


# Calculate the total elapsed time for creating and submitting all job scripts
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

# Print the total time taken
echo "Job scripts created and submitted in $elapsed_time seconds."

#################################################################  Longread Island Construction    ##############################################################################


python3 $HOME/Maptcha/src/CreateUnmappedUnusedLR.py $HOME/Maptcha/Output/FastaFilesBatch_8192/ $HOME/Maptcha/Output/contExp.fasta $HOME/Maptcha/TestInput/CoxiellaBurnetii_longreads.fa $HOME/Maptcha/Output/unused_longreads.fasta


## Define compute options
#PBS -l nodes=3:amd:ppn=20
#PBS -l mem=120gb
#PBS -l walltime=06:00:00


cd $HOME/Maptcha/Hifiasm/
chmod +x $HOME/Maptcha/Hifiasm/hifiasm

# Start the timer
start_time=$(date +%s)

mkdir $HOME/Maptcha/Output/Phase2/

$HOME/Maptcha/Hifiasm/hifiasm -o $HOME/Maptcha/Output/Phase2/Only_UnmappedUnusedLongreads.asm -t 64 $HOME/Maptcha/Output/unused_longreads.fasta

awk '/^S/{print ">"$2;print $3}' $HOME/Maptcha/Output/Phase2/Only_UnmappedUnusedLongreads.asm.bp.p_ctg.gfa > $HOME/Maptcha/Output/Phase2/Only_UnmappedUnusedLongreads.asm.bp.p_ctg.gfa.fa


# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

echo "Longread Island Construction done! "
#################################################################  Link Scaffolds With Bridges    ##############################################################################


python3 $HOME/Maptcha/src/merge.py $HOME/Maptcha/Output/Phase2/Only_UnmappedUnusedLongreads.asm.bp.p_ctg.gfa.fa $HOME/Maptcha/Output/contExp.fasta $HOME/Maptcha/Output/Phase1_2_partialScaff.fa


## Define compute options
#PBS -l nodes=1:amd:ppn=20
#PBS -l mem=120gb
#PBS -l walltime=06:00:00


cd $HOME/Maptcha/Hifiasm/
chmod +x $HOME/Maptcha/Hifiasm/hifiasm

# Start the timer
start_time=$(date +%s)

mkdir $HOME/Maptcha/Output/Final/

$HOME/Maptcha/Hifiasm/hifiasm -o $HOME/Maptcha/Output/Final/finalAssembly.asm -t 64 -n 1 "$phase1_2_output" "$unusedlongreads"

awk '/^S/{print ">"$2;print $3}' $HOME/Maptcha/Output/Final/finalAssembly.asm.bp.p_ctg.gfa > $HOME/Maptcha/Output/Final/finalAssembly.fa


# Calculate the elapsed time
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
echo "Link Scaffolds With Bridges done! "


