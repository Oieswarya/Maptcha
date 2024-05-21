// ***********************************************************************
//
// PaKman: Algorithm for generating genomic contigs on distributed-memory machines
// 
// Priyanka Ghosh (Pacific Northwest National Laboratory)
// Sriram Krishnamoorthy (Pacific Northwest National Laboratory)
// Ananth Kalyanaraman (Washington State University)
//               
//
// ***********************************************************************
//
//       Copyright (2020) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "JEM.h"
#include "timers.h"
#include <iostream>
#include <fstream>
#include <sstream>
//#include "serial.h"

using std::cout; using std::cerr;
using std::endl; using std::string;
using std::ifstream; using std::ostringstream;

long int MAX_KMER_COUNT=0;
int rank, size;
int coverage=0;
std::string inputFileName;
std::string queryFileName;
std::string primeFileName;
std::string AFileName;
std::string BFileName;
int read_length=0;

int node_threashold=0;

int num_batch_transfers=0;


//std::vector<lmer_t> lmer_frequency(LMER_SIZE,0);
//std::vector<lmer_t> global_lmer_frequency(LMER_SIZE,0);

/* Vector containing k-mers allocated to each process */
std::vector<KmerPairs> kmer_proc_buf;
std::vector<std::vector<kmer_t>> kmer_sets;
//std::vector<std::unordered_map<kmer_t, std::vector<int>> > Tl;

//kmer_t *hasha, *hashb;
uint64_t hasha = 68111;
uint64_t hashb = 105929;




void parseCommandLine(const int argc, char * const argv[]);




std::string readFileIntoString(const std::string& path) {
    ifstream input_file(path);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << path << "'" << endl;
        exit(EXIT_FAILURE);
    }
    return string((std::istreambuf_iterator<char>(input_file)), std::istreambuf_iterator<char>());
}

int get_file_size(std::string filename) // path to file
{
    FILE *p_file = NULL;
    p_file = fopen(filename.c_str(),"rb");
    fseek(p_file,0,SEEK_END);
    int size = ftell(p_file);
    fclose(p_file);
    return size;
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    parseCommandLine(argc, argv);
     

    //std::string contigFileName = inputFileName + std::to_string(rank) + ".fa";
    //std::cout<<std::to_string(rank)<<" - "<<inputName<<"\n";
    double time_l1 = MPI_Wtime ();
    
    //input_read_data cdata = perform_input_reading(rank, size, inputFileName, 108539);
   // input_read_data cdata = perform_input_reading(rank, size, inputFileName, 25980); //34680 //25980 EC //27586 PA
    input_read_data cdata = perform_input_reading(rank, size, inputFileName, 127978);
    //input_read_data rdata = perform_input_readingC(rank, size, inputName, read_length);
    
    //string file_contents;
    //size_t file_size;
    
    //std::cout<<"Print statement\n";
    //std::cout<<contigFileName<<"\n";
        //contigdata = perform_input_readingC(rank, size, contigFileName, read_length);
    //file_contents = readFileIntoString(contigFileName);
    //file_size = get_file_size(contigFileName);
    //std::cout<<"Print statement\n";
        
        //file_contents = readFileIntoString(filename);
    
    //read_array();
    //perform_kmer_counting (file_contents, file_size);
    int M;
    int total_subjects;
    
    
    input_read_data rdata = perform_input_reading(rank, size, queryFileName, read_length);
    double time_l3 = MPI_Wtime ();
    generate_set_of_subjects (cdata.read_data, cdata.read_data_size, cdata.start_index, rdata.read_data, rdata.read_data_size, rdata.start_index, &M, &total_subjects);
    //int total_number_of_subs_in_p = kmer_sets.size();
    //kmer_t** Hash_table = new kmer_t*[total_subjects];
    
    /*int total_hash_functions = 150;
    for (int i = 0; i < total_subjects; i++) {
        Hash_table[i] = new kmer_t[total_hash_functions];
    }
    genereate_hash_table(M, total_subjects, Hash_table);*/
    //double time_l11 = MPI_Wtime ();
    
    double time_l2 = MPI_Wtime ();
    double f_time = time_l2 - time_l3;
    //printf ("%d Average time for out across all procs (secs): %f \n", rank, f_time);
    //free(cdata.read_data);
    //free(rdata.read_data);
    /*
    if(rank == 0)
    {
        const std::string s0("/home/trahman/GitRepo/pakman-pakman_latest/Hash_");
        char proc_id[3];
        char output_file_name[25];

        sprintf(proc_id, "%d", rank); 
    //strcpy(output_file_name,"wired_mn_out_");
        strcpy(output_file_name, s0.c_str());
        strcpy(&output_file_name[strlen(output_file_name)],proc_id);
        strcpy(&output_file_name[strlen(output_file_name)],".log");
        FILE *f = fopen(output_file_name, "w");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        for(int h = 0; h<total_subjects; h++)
        {
            for(int k = 0; k<200; k++)
            {
                //fprintf(f,"%ld ", Hash_table[h][k]); 
            }
        //    fprintf(f, "\n");
        }
    }
    
    */
    //double time_l3 = MPI_Wtime ();
    //input_read_data rdata = perform_input_reading(rank, size, queryFileName, read_length);


    //perform_kmer_counting (rdata.read_data, rdata.read_data_size);
   // generate_set_of_queries (rdata.read_data, rdata.read_data_size, rdata.start_index, total_subjects, M, Hash_table);
    
    double time_l4 = MPI_Wtime ();
    double s_time = time_l4 - time_l3;
    //printf ("%d Average time for l across all procs (secs): %f \n", rank, s_time);

    MPI_Finalize();
    return 0;
}

void parseCommandLine(const int argc, char * const argv[])
{
  int ret;

  while ((ret = getopt(argc, argv, "s:q:p:a:b:l:n:")) != -1) {
    switch (ret) {
    case 's':
       inputFileName.assign(optarg);
       //std::cout << inputFileName << std::endl;
       break;
    case 'q':
       queryFileName.assign(optarg);
       //std::cout << inputFileName << std::endl;
       break;
    case 'p':
       primeFileName.assign(optarg);
       //std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'a':
       AFileName.assign(optarg);
       //std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'b':
       BFileName.assign(optarg);
       //std::cout << MAX_KMER_COUNT << std::endl;
       break;
    case 'l':
       read_length = atoi(optarg);
       //std::cout << read_length << std::endl;
       break;
    
    case 'n':
       node_threashold = atoi(optarg);
       //std::cout << node_threashold << std::endl;
       break;
    default:
       assert(0 && "Should not reach here!!");
       break;
    }
  }

  if (rank==0 && inputFileName.empty()) {
      std::cerr << "Must specify an input FASTA file name with -f" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

//  if (rank==0 && !MAX_KMER_COUNT) {
//      std::cerr << "Must specify a batch size for k-mer counting with -b" << std::endl;
//      MPI_Abort(MPI_COMM_WORLD, -99);
//  }

  

  if (rank==0 && !read_length) {
      std::cerr << "Must specify read_length with -r" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -99);
  }

  

  


} // parseCommandLine


