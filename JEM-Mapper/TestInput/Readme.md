# Generate contigs using Minia assembler:

To generates contigs from short reads we have used minia short reads assembler (https://github.com/GATB/minia).
Please find their complete manual: https://github.com/GATB/minia/raw/master/doc/manual.pdf

The basic usage is as follows where input_file is the input short read fasta file:

./minia -in [input file] -kmer-size [kmer size] -out [prefix]

To generate contigs from Minia assembler, we have used kmer-size as 32.


# Generate long reads using Sim-it HiFi read simulatr:

We have generated HiFi long reads using Sim-it simulator (https://github.com/ndierckx/Sim-it). 

To generate 10x coverage HiFi long reads, we made the following changes to the configuration file.
* Sequencing depth         = 10
* Median length            = 10000
* Length range             = 10000-15000
* Accuracy                 = 99.9%

Since, we wanted 10x coverage long reads, we set sequncing depth as 10, we set median length as 10,000 and length variance was expected  from 10,000 to 15,000.
And since, we expected 99.9% accuracy, we set Accuracy as 99.9%.
