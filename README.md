# Maptcha
Maptcha addresses the hybrid scaffolding problem. We have three major steps: (1) Mapping the contigs to long reads, (2) Constructing the contig-contig graph from the mapping information, and (3) Ordering the contigs.

For the first step, mapping, we use JEM-Mapper: https://github.com/TazinRahman1105050/JEM-Mapper
Please use the link to install and map your contigs to long reads.

For the graph construction, we use a python script: graphconstr.py
For ordering the contigs there are three approaches we used---(i)b-way graph matching (https://github.com/ECP-ExaGraph/bMatching), (ii) wiring heuristic (Wiring.py), and (iii) vertex reordering (CreateRCM.py).

The combined.py script contains all the three python scripts.

toy.log is for testing. Contains the graph. Can be used as input for wiring heusristic and vertex reordering approach.

toyForBway.mtx is test input for bway graph matching.
