# Maptcha

Maptcha addresses the hybrid scaffolding problem. We have three major phases: 
1. Contig Expansion: Initially, the algorithm extends contigs using long reads that align with the ends of these contigs. This phase also involves detecting and connecting successive pairs of contigs using direct long read links, resulting in the generation of partial scaffolds.

2. Long Read Island Construction: Not all long reads contribute to the initial partial scaffolds, especially those residing in the gap regions between successive scaffolds in the target genome. In this phase, the algorithm identifies long reads that do not map to any first-generation partial scaffolds. These reads are utilized to construct partial scaffolds corresponding to the "island" regions of long reads, forming the second generation of partial scaffolds.

3. Link Scaffolds with Bridges: In the final phase, the algorithm aims to bridge the first and second generation scaffolds using long reads that serve as bridges between them. This crucial step produces the final set of scaffolds, providing a comprehensive assembly of the genome.



## Installation Instructions

Maptcha utilizes the following tools:

- **JEM-Mapper**: [JEM-Mapper GitHub Repository](https://github.com/TazinRahman1105050/JEM-Mapper)
- **Hifiasm**: [Hifiasm GitHub Repository](https://github.com/chhylp123/hifiasm)

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

### Usage
Run the maptcha.sh script from the root directory:

```bash
./maptcha.sh -c path/to/contigs.fa -lr path/to/longreads.fa [options]

Options:
-c, --contigs Path to the contigs input file
-lr, --longreads Path to the long reads input file
-t, --threads Number of threads to use (default: 1)
-n, --nodes Number of nodes to use (default: 2)
-p, --processes Number of processes per node (default: 2)
-h, --help Show this help message
Notes:
Makefile: Manages setup, compilation, permissions, cleaning, and dependency installation.
requirements.txt: Lists Python dependencies required for your tool (biopython, networkx, tqdm).
README.md: Provides instructions on how to clone, set up dependencies, compile, and use your tool.


### For a quick test, you can use the provided setup where the necessary binary of the tools are pre-installed for a test input. Navigate within the Maptcha repository and run the `maptcha.sh` script. 

```bash
~/Maptcha/src/maptcha.sh ~/Maptcha/TestInput/minia_Coxiellaburnetii_contigs.fa ~/Maptcha/TestInput/CoxiellaBurnetii_longreads.fa
```
Ensure that you have the appropriate permissions to execute the job script.

The final scaffolds will be located here: `~/Maptcha/Output/Final/finalAssembly.fa`, within the Output folder of the Maptcha directory.

**Requirements:**
You will still need the following modules:
- C++14 (or greater) compliant compiler
- MPI library (preferably MPI-3 compatible)
- Python 3 (or greater) compliant compiler:

For more detailed usage and configuration options, please refer to the documentation within each tool's repository:

- [JEM-Mapper Documentation](https://github.com/TazinRahman1105050/JEM-Mapper)
- [Hifiasm Documentation](https://github.com/chhylp123/hifiasm)

If you encounter any issues, please feel free to open an issue on the GitHub repository.
