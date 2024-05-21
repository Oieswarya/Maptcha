# Maptcha

Maptcha addresses the hybrid scaffolding problem. We have three major phases: 
1. Contig Expansion: Initially, the algorithm extends contigs using long reads that align with the ends of these contigs. This phase also involves detecting and connecting successive pairs of contigs using direct long read links, resulting in the generation of partial scaffolds.

2. Long Read Island Construction: Not all long reads contribute to the initial partial scaffolds, especially those residing in the gap regions between successive scaffolds in the target genome. In this phase, the algorithm identifies long reads that do not map to any first-generation partial scaffolds. These reads are utilized to construct partial scaffolds corresponding to the "island" regions of long reads, forming the second generation of partial scaffolds.

3. Link Scaffolds with Bridges: In the final phase, the algorithm aims to bridge the first and second generation scaffolds using long reads that serve as bridges between them. This crucial step produces the final set of scaffolds, providing a comprehensive assembly of the genome.

################################################################################
**Requirements:**
- C++14 (or greater) compliant compiler
- MPI library (preferably MPI-3 compatible)
- Python 3 (or greater) compliant compiler
  
################################################################################

# Maptcha

## Installation Instructions

To map long reads to contigs, Maptcha utilizes the following tools:

- **JEM-Mapper**: [JEM-Mapper GitHub Repository](https://github.com/TazinRahman1105050/JEM-Mapper)
- **Hifiasm**: [Hifiasm GitHub Repository](https://github.com/chhylp123/hifiasm)

### Step-by-Step Guide

1. **Install JEM-Mapper and Hifiasm**

   Both JEM-Mapper and Hifiasm have their own dependencies and installation procedures. Please refer to their respective repositories for detailed installation instructions. We recommend installing these tools in the same directory as Maptcha to ensure smooth integration.

2. **Clone the Maptcha Repository**

   ```bash
   git clone https://github.com/Oieswarya/Maptcha.git

### For a quick test, you can use the provided setup with a test input where the necessary tools are pre-installed. Navigate  within the Maptcha repository and run the `maptcha.sh` script:


```bash
chmod +x maptcha/src/maptcha.sh
./maptcha.sh
```

Ensure that you have the appropriate permissions to execute the job script.

For more detailed usage and configuration options, please refer to the documentation within each tool's repository:

- [JEM-Mapper Documentation](https://github.com/TazinRahman1105050/JEM-Mapper)
- [Hifiasm Documentation](https://github.com/chhylp123/hifiasm)

By following these steps, you should be able to set up and run Maptcha seamlessly. If you encounter any issues, please feel free to open an issue on the GitHub repository.
