# Maptcha
Maptcha addresses the hybrid scaffolding problem. We have three major phases: 

1. Contig Expansion: Initially, the algorithm extends contigs using long reads that align with the ends of these contigs. This phase also involves detecting and connecting successive pairs of contigs using direct long read links, resulting in the generation of partial scaffolds.

2. Long Read Island Construction: Not all long reads contribute to the initial partial scaffolds, especially those residing in the gap regions between successive scaffolds in the target genome. In this phase, the algorithm identifies long reads that do not map to any first-generation partial scaffolds. These reads are utilized to construct partial scaffolds corresponding to the "island" regions of long reads, forming the second generation of partial scaffolds.

3. Link Scaffolds with Bridges: In the final phase, the algorithm aims to bridge the first and second generation scaffolds using long reads that serve as bridges between them. This crucial step produces the final set of scaffolds, providing a comprehensive assembly of the genome.

In the the first phase, for the purpose of mapping the long reads to contigs, we use JEM-Mapper:(https://github.com/TazinRahman1105050/JEM-Mapper)]
Please use the link to install and map your contigs to long reads.

For the graph construction, we use a python script: graphconstr.py
For ordering the contigs there are three approaches we used---(i)b-way graph matching (https://github.com/ECP-ExaGraph/bMatching), (ii) wiring heuristic (Wiring.py), and (iii) vertex reordering (CreateRCM.py).

The combined.py script contains all the three python scripts.

Please read the ToyInputData/ReadmeForInputs.md to use the sample inputs efficiently.
