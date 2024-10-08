# Maptcha

Maptcha addresses the hybrid genome scaffolding problem, which involves combining contigs and long reads to create a more complete and accurate genome assembly. Hybrid scaffolding aims to leverage the high accuracy of contigs with the long-range information provided by long reads.

Maptcha's inputs are a FASTA file of contigs and a FASTA file of long reads, and the output is a FASTA file of scaffolds generated from them.

We have three major phases: 
1. Contig Expansion: Initially, the algorithm extends contigs using long reads that align with the ends of these contigs. This phase also involves detecting and connecting successive pairs of contigs using direct long read links, resulting in the generation of partial scaffolds.

2. Long Read Island Construction: Not all long reads contribute to the initial partial scaffolds, especially those residing in the gap regions between successive scaffolds in the target genome. In this phase, the algorithm identifies long reads that do not map to any first-generation partial scaffolds. These reads are utilized to construct partial scaffolds corresponding to the "island" regions of long reads, forming the second generation of partial scaffolds.

3. Link Scaffolds with Bridges: In the final phase, the algorithm aims to bridge the first and second generation scaffolds using long reads that serve as bridges between them. This crucial step produces the final set of scaffolds, providing a comprehensive assembly of the genome.

## Citation
If you use Maptcha in your research, please cite:

Bhowmik, O., Rahman, T. & Kalyanaraman, A. Maptcha: an efficient parallel workflow for hybrid genome scaffolding. BMC Bioinformatics 25, 263 (2024). https://doi.org/10.1186/s12859-024-05878-4

## Installation Instructions

**Requirements:**
Maptcha has the following dependencies:
- C++14 (or greater) compliant compiler         
- MPI library (MPI-3 compatible)     
- Python 3 (or greater) compliant compiler      


### Step-by-Step Guide

1. **Clone the Maptcha Repository:**

   ```bash
   git clone https://github.com/Oieswarya/Maptcha.git
   cd Maptcha
   
2. **Install necessary dependencies:**

   ```bash
   make install-dependencies

3. **Compile the source files and setup directories:**

   ```bash
   make all

3. **Check if Maptcha is properly installed:**

   ```bash
   ./maptcha.sh -h

### Usage
Run the maptcha.sh script from the root directory:

```bash
./maptcha.sh -c path/to/contigs.fa -lr path/to/longreads.fa [options]

-c, --contigs      Path to the contigs input file
-lr,--longreads    Path to the long reads input file
Options:
-o, --output       Output directory (default: $HOME/Maptcha/Output/)
-t, --threads      Number of threads to use (default: 16)
-n, --nodes        Number of nodes to use (default: 2)
-p, --processes    Number of processes per node (default: 2)
-h, --help         Show this help message
```

Note:
This code has been tested on high-performance cluster (HPC) systems with MPI and OpenMP compatibility and has been tested for both PBS and SLURM job scheduling systems.


### For a quick test, you can use the provided test input. Navigate within the Maptcha repository and run the `maptcha.sh` script. 

```bash
~/Maptcha/maptcha.sh -c ~/Maptcha/TestInput/minia_Coxiellaburnetii_contigs.fa -lr ~/Maptcha/TestInput/CoxiellaBurnetii_longreads.fa
```

The final scaffolds will be located here: `~/Maptcha/Output/Final/finalAssembly.fa`, within the Output folder of the Maptcha directory.


**Tips:**
1. On some clusters, you may need to load specific modules before installing dependencies and and then also while running Maptcha.
2. Ensure that you have the appropriate permissions to execute the job script.

Maptcha utilizes the following tools:

- **JEM-Mapper**: [JEM-Mapper GitHub Repository](https://github.com/TazinRahman1105050/JEM-Mapper)
- **Hifiasm**: [Hifiasm GitHub Repository](https://github.com/chhylp123/hifiasm)

For more detailed usage and configuration options, please refer to the documentation within each tool's repository:

- [JEM-Mapper Documentation](https://github.com/TazinRahman1105050/JEM-Mapper)
- [Hifiasm Documentation](https://github.com/chhylp123/hifiasm)

If you encounter any issues, please feel free to open an issue on the GitHub repository.
