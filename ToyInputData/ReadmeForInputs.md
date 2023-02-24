We have added a toy input contigs and LR fasta file. The first step is to map the contigs to the long redas. For this step, we have used JEM-Mapper.
The output of JEM-Mapper is the contigs mapped to the long reads.

The next step is the graph construction. The code is written in the graphconstr.py. 
The graph generated from this step and the sample fasta files will look like this:

1 2 2

2 3 1

3 4 2

4 5 2

This is the input for the scaffolding techniques. 

We have also added the log files which are essentially the output after the graph construction step, in case, anyone wants to use it directly for the scaffolding steps.
