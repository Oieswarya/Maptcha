#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <omp.h>

using namespace std;

struct ContigPair {
    int contig_id1;
    int contig_id2;
    int longread_id;
};

bool compareContigPair(const ContigPair& cp1, const ContigPair& cp2) {
    if (cp1.contig_id1 < cp2.contig_id1) {
        return true;
    } else if (cp1.contig_id1 > cp2.contig_id1) {
        return false;
    } else {
        return cp1.contig_id2 < cp2.contig_id2;
    }
}

int main(int argc, char** argv) {
    double t1,t2,t3,t4,t5;
    t1=omp_get_wtime();
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input_file output_file\n";
        exit(1);
    }

    string input_filename = argv[1];
    string output_filename = argv[2];

    ifstream input_file(input_filename);
    if (!input_file) {
        cerr << "Error: Could not open input file " << input_filename << endl;
        exit(1);
    }

    ofstream output_file(output_filename);
    if (!output_file) {
        cerr << "Error: Could not open output file " << output_filename << endl;
        exit(1);
    }

    vector<pair<int, int>> contig_longreads;
    int contig_id, longread_id;
    while (input_file >> contig_id >> longread_id) {
        contig_longreads.emplace_back(contig_id, longread_id);
    }
    input_file.close();

    sort(contig_longreads.begin(), contig_longreads.end());

    unordered_map<int, vector<int>> longreads_contigs;
    for (const auto& cl : contig_longreads) {
        longreads_contigs[cl.second].push_back(cl.first);
    }

    vector<ContigPair> contig_pairs;
    #pragma omp parallel for
    for (int i = 0; i < contig_longreads.size(); i++) {
        const auto& cl = contig_longreads[i];
        int contig_id1 = cl.first;
        int longread_id = cl.second;

        auto& contigs_with_same_longread = longreads_contigs[longread_id];
        int start_index = lower_bound(contigs_with_same_longread.begin(), contigs_with_same_longread.end(), contig_id1+1) - contigs_with_same_longread.begin();
        for (int j = start_index; j < contigs_with_same_longread.size(); j++) {
            int contig_id2 = contigs_with_same_longread[j];
            #pragma omp critical
            {
                contig_pairs.push_back({contig_id1, contig_id2, longread_id});
            }
        }
        
    }
    
    sort(contig_pairs.begin(), contig_pairs.end(), compareContigPair);
    
    for (const auto& cp : contig_pairs) {
        output_file << cp.contig_id1 << " " << cp.contig_id2 << " " << cp.longread_id << endl;
    }
    
    output_file.close();
    
    t2=omp_get_wtime();
    cout << "Total Time: " << t2-t1 << endl;

    return 0;
}
