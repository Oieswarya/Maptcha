#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <omp.h>

using namespace std;

int main(int argc, char** argv) {
    double t1, t2;
    t1 = omp_get_wtime();

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

    // Sorting contig_longreads by long read IDs
    sort(contig_longreads.begin(), contig_longreads.end(),
         [](const pair<int, int>& a, const pair<int, int>& b) {
             return a.second < b.second;
         });

    unordered_map<int, unordered_map<int, int>> contig_pairs_count;
#pragma omp parallel for
    for (int i = 0; i < contig_longreads.size(); i++) {
        int contig_id1 = contig_longreads[i].first;
        int longread_id = contig_longreads[i].second;

        for (int j = i + 1; j < contig_longreads.size() && contig_longreads[j].second == longread_id; j++) {
            int contig_id2 = contig_longreads[j].first;
            contig_pairs_count[contig_id1][contig_id2]++;
        }
    }

    // Output generation and writing
    int max_contig_id = 0;
    for (const auto& cl : contig_longreads) {
        max_contig_id = max(max_contig_id, cl.first);
    }

    int num_pairs = 0;
    for (const auto& cp : contig_pairs_count) {
        num_pairs += cp.second.size();
    }

    output_file << "%%MatrixMarket matrix coordinate real symmetric" << endl;
    output_file << max_contig_id << " " << max_contig_id << " " << num_pairs << endl;

    for (int contig_id1 = 0; contig_id1 <= max_contig_id; contig_id1++) {
        const auto it = contig_pairs_count.find(contig_id1);
        if (it == contig_pairs_count.end()) {
            continue;
        }

        const auto& pairs_count = it->second;
        for (const auto& pc : pairs_count) {
            int contig_id2 = pc.first;
            int count = pc.second;
            if (count > 0) {
                output_file << contig_id1 << " " << contig_id2 << " " << count << endl;
            }
        }
    }

    output_file.close();

    t2 = omp_get_wtime();
    cout << "Total Time: " << t2 - t1 << endl;

    return 0;
}
