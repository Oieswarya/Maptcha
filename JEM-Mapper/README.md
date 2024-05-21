# Dependencies:
JEM-Mapper has the following dependencies:

* MPI library (preferably MPI-3 compatible)
* C++14 (or greater) compliant compiler

# Build:
make ksize=$KMER_SIZE

For example:
make ksize = 15

# Execute:
mpiexec -np $number_of_procs $BINARY -s {Contig_Fasta_File} -q {Long_Read_Fasta_File} -a {A_int_Values_File} -b {B_int_Values_File} -p {Prime_int_Values_File} -r $read_segment length -n $NO_OF_TRIALS

For example:

mpiexec -np 4 ./jem -s ~/Ecoli_reads_100x_contigs.fasta -q ~/Ecoli_reads_10x_long_reads.fasta -a ~/A.txt -b ~/B.txt -p ~/Prime.txt -r 1000 -n 30

Notes:
* This code has been tested on high-performance computing cluster (HPC) with MPI compatibility. For the system we used we had to set the number of processes in the given way. Please change the parameters accordingly.

# Input arguments 
* -s: input contigs fasta file
* -q: input long reads fasta file
* -a: For diffent trials we have used a linear congruential hash funtion of the form: (Ax+B)%P, whwere A and B are integers and P is a prime and x is a kmer we want to hash. So, using this parameter, we provide the input values for different A values
* -b: input B values
* -p: input prime numbers
* -r: read segment length
* -n: number of trials
