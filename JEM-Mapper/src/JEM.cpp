
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
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
#include <fstream>
using std::ifstream; using std::ostringstream;

extern int rank, size;
extern int num_threads;
extern int node_threashold;
extern std::string primeFileName;
extern std::string AFileName;
extern std::string BFileName;

int Ax[200];
int Bx[200];
int Px[200];

extern std::vector<std::vector<kmer_t>> kmer_sets;

//extern std::vector<std::unordered_map<kmer_t, std::vector<int>> > Tl(10);


void read_array()
{
    std::ifstream input(primeFileName);

    for (int i = 0; i < 200; i++) {
        input >> Px[i];
        
        //std::cout<< A[i]<<std::endl;
    }
    std::ifstream input1(AFileName);

    for (int i = 0; i < 200; i++) {
        input1 >> Ax[i];
        //std::cout<< A[i]<<std::endl;
    }
    std::ifstream input2(BFileName);

    for (int i = 0; i < 200; i++) {
        input2 >> Bx[i];
        //std::cout<< A[i]<<std::endl;
    }
}



void generate_modified_set(int M, std::vector<std::vector<kmer_t>> &previous_sets)
{
    kmer_t start_kmer_for_modified_set = (1UL<<(2*16)); //C followed by AAA...s
    //printf("%d () %ld\n",rank, start_kmer_for_modified_set);
    //int update = 0;
    for(int i = 0; i<previous_sets.size(); i++)
    {
        int Modified = M - previous_sets[i].size();
        for(int j = 0; j< Modified; j++)
        {
            previous_sets[i].push_back(start_kmer_for_modified_set+j);
            
        }
      /*  if(Modified > 0)
        {
            update++;
        }*/
    }
}

void generate_modified_set_queries(int M, std::vector<std::vector<kmer_t>> &previous_sets)
{
   // std::cout<<"\n";
    kmer_t start_kmer_for_modified_set = (1UL<<(2*16)); //C followed by AAA...s
    //printf("%d () %ld\n",rank, start_kmer_for_modified_set);
    int update = 0;
    for(int i = 0; i<previous_sets.size(); i++)
    {
        if(previous_sets.size() > 0)
        {
        int Modified = M - previous_sets[i].size();
        for(int j = 0; j< Modified; j++)
        {
            previous_sets[i].push_back(start_kmer_for_modified_set+j+M);
            
        }
        if(Modified > 0)
        {
            update++;
        }
        }
    }
    //printf("%d () %ld\n",rank, update);
}
/*
void recalculate_min_lmer (kmer_t kmer_in, lmer_t *m_lmer, lmer_t *m_lmer_freq, int *m_pos)
{
    lmer_t min_lmer=0, tmp_lmer=0;
    lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, k=0;

    for (k=0; ((10-1) - k) >= (16-1); k++) {
        lmer_t lmer_out=0;
        for(int j=k; j<16+k; j++) {
            lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
        }

        tmp_lmer = lmer_out;
        //tmp_lmer_freq = global_lmer_frequency[tmp_lmer];

        if (k == 0) {
            min_lmer = tmp_lmer;
            min_lmer_freq = tmp_lmer_freq;
            min_pos = 0;
        }
        else {
           if (tmp_lmer < min_lmer) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
               min_pos = k;
           }
        }
    }
    assert (k == (10-16+1));

    *m_lmer = min_lmer;
    *m_lmer_freq = min_lmer_freq;
    *m_pos = min_pos;
}
*/
void recalculate_min_kmer (std::string ptr, kmer_t *m_kmer, int *fact, int *pos)
{
    kmer_t min_lmer=3074457345618258602, tmp_lmer=0;
    //lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    //int min_pos=0, k=0;

    kmer_t kmer=0;
    //rev_kmer=0;
    int i;
    //std::string s;
    int tracker = 0;
    //int contig_len = 0;
    int tt = 0; //tracking N
    //std::cout << ptr<<"\n";
    for(int i=0; i<16-1; i++) {
            //N Cheking
        if(tt > 0)
        {
            tt--;
        }
        if (ptr[i] == 'N' || ptr[i] == 'Y' || ptr[i] == 'S' || ptr[i] == 'R' || ptr[i] == 'I')
        {
            tt = 16;
        }
            //N Checking
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
        kmer = kmer_shift(kmer, char_to_el(ptr[i]));
            /*
            if(rank == 0)
            {
                std::cout<<ptr[p]<<"\n";
            }*/
        //    s.push_back(convert_to_char(ptr[p]));
        //    p++;
        //    contig_len++;
    }
    int start = 0;
    for(int i=16-1; i < 100; i++) {
    //while(p<length && !isspace(ptr[p])) {
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //N Cheking
            if(tt > 0)
            {
                tt--;
            }
            if (ptr[i] == 'N' || ptr[i] == 'Y' || ptr[i] == 'S' || ptr[i] == 'R' || ptr[i] == 'I')
            {
                tt = 16;
            }
            //N Checking
            
            kmer = kmer_shift(kmer, char_to_el(ptr[i]));
            /*if (i == 16 -1)
            {
                
                printf("+-%ld-+\n", kmer); 
            }*/
            //recalculate_min_lmer(kmer, &min_kmer, &min_lmer_freq, &min_pos);
            /*if(rank == 0 && chkr == 0)
            {
                //printf("rank %d , %ld\n", rank, kmer);
                chkr += 1;
            }*/
            //s.push_back(convert_to_char(ptr[p]));
            //s.push_back(ptr[p]);
            //p++;
            //contig_len++;
            //rev_set.push_back(kmer);
            /*if (tt > 0)
            {
                //rev_set_tracker.push_back(0);
                //tracker = 0; 
            }*/
            if(tt <= 0)
            {
                tracker = 1;
                //rev_set_tracker.push_back(1);
                if (min_lmer > kmer)
                {
                    min_lmer = kmer;
                    start = i - 15;
                }
            }
            //start += 1;
        } 
    *m_kmer = min_lmer;
    //*m_lmer_freq = min_lmer_freq;
    *fact = tracker;
    *pos = start;
}

/*
void get_hash_value(std::vector<kmer_t> previous_sets, int *A1, int *A2, int *A3, int *A4, int *A5)
{
    int min1=2147483647, min2=2147483647, min3=2147483647, min4=2147483647, min5=2147483647;
    for(int i = 0; i<previous_sets.size(); i++)
    {
        int val = previous_sets[i]%27644437;
        if(val<min1)
        {
            min1 = val;
        }
        int val1 = previous_sets[i]%999331;
        if(val1<min2)
        {
            min2 = val1;
        }
        int val2 = previous_sets[i]%319993;
        if(val2<min3)
        {
            min3 = val2;
        }
        int val3 = previous_sets[i]%116109;
        if(val3<min4)
        {
            min4 = val3;
        }
        int val4 = previous_sets[i]%933199;
        if(val4<min5)
        {
            min5 = val4;
        }
    }
    *A1 = min1;
    *A2 = min2;
    *A3 = min3;
    *A4 = min4;
    *A5 = min5;
}
*/


void prefix_val(kmer_t **A1, int M)
{
    //for(int j = 0; j<kmer_sets.size(); j++)
    //{
    //int min1=2147483647, min2=2147483647, min3=2147483647, min4=2147483647, min5=2147483647;
        kmer_t Max_kmer_val = 3074457345618258602;
        std::vector<kmer_t> minL (100,Max_kmer_val);
        std::vector<kmer_t> corr_kmer (100,3074457345618258602);
        kmer_t start_kmer_for_modified_set = (1UL<<(2*KMER_LENGTH));
        for(int i = 0; i<M; i++)
        {
            kmer_t new_k = start_kmer_for_modified_set+i;
            for (int k = 0; k < 100; k++)
            {
               // if(j == 0 && rank == 0 && i == 0)
                //{   
          //      printf("%d %ld %ld %ld %ld\n",rank, kmer_sets[j][i], Ax[k], Bx[k], Px[k]);
          //      printf("%d %ld\n",rank, (Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
            //}
                //kmer_t val = (((Ax[k]%Px[k])*(kmer_sets[j][i]%Px[k]))%Px[k] + (Bx[k]%Px[k]))%Px[k];
                 kmer_t val = (Ax[k]*new_k + Bx[k]) % Px[k];
            //if (((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]) < minL[k])
                if (val <minL[k])
                {
                //minL[k] = ((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
                    minL[k] = val;
                    corr_kmer[k] = new_k;
                }
                A1[k][i] = corr_kmer[k];
            }
        
        }
        
        corr_kmer.clear();
        corr_kmer.shrink_to_fit();
        minL.clear();
        minL.shrink_to_fit();
    //}
    //kmer_sets.clear();
    //kmer_sets.shrink_to_fit();
}

void get_hash_value(kmer_t **A1, int M, kmer_t **Prefix)
{
    for(int j = 0; j<kmer_sets.size(); j++)
    {
    //int min1=2147483647, min2=2147483647, min3=2147483647, min4=2147483647, min5=2147483647;
        kmer_t Max_kmer_val = 3074457345618258602;
        std::vector<kmer_t> minL (100,Max_kmer_val);
        std::vector<kmer_t> corr_kmer (100,3074457345618258602);
    
        for(int i = 0; i<kmer_sets[j].size(); i++)
        {
            for (int k = 0; k < 100; k++)
            {
               // if(j == 0 && rank == 0 && i == 0)
                //{   
          //      printf("%d %ld %ld %ld %ld\n",rank, kmer_sets[j][i], Ax[k], Bx[k], Px[k]);
          //      printf("%d %ld\n",rank, (Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
            //}
                //in
                //kmer_t val = (((Ax[k]%Px[k])*(kmer_sets[j][i]%Px[k]))%Px[k] + (Bx[k]%Px[k]))%Px[k];
                 kmer_t val = (Ax[k]*kmer_sets[j][i] + Bx[k]) % Px[k];
            //if (((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]) < minL[k])
                if (val <minL[k])
                {
                //minL[k] = ((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
                    minL[k] = val;
                    corr_kmer[k] = kmer_sets[j][i];
                }
            }
        
        }
        for (int k = 0; k < 100; k++)
        {
            kmer_t prefix = Prefix[k][M-kmer_sets[j].size()-1];
            kmer_t p_val = (Ax[k]*prefix + Bx[k]) % Px[k];
            if (minL[k] <= p_val)
            {
                A1[j][k] = corr_kmer[k];
            }
            else
            {
                A1[j][k] = prefix;
            }
            //if(j == 0 && rank == 0)
            //{
        //printf("%d %ld\n", rank, A1[j][k]);
            //}
    //A1[j][1] = min2;
    //A1[j][2] = min3;
    //A1[j][3] = min4;
    //A1[j][4] = min5;
        }
        corr_kmer.clear();
        corr_kmer.shrink_to_fit();
        minL.clear();
        minL.shrink_to_fit();
    }
    kmer_sets.clear();
    kmer_sets.shrink_to_fit();
}

void get_hash_value_queires1(std::vector<std::vector<kmer_t>> &modified_sets, kmer_t **A1)
{
    //double oq, tq, oq1, oq2 = 0.0;
    kmer_t val;
    //std::cout<<oq<< " "<<tq<<"\n";
    for(int j = 0; j<modified_sets.size(); j++)
    {
    //int min1=2147483647, min2=2147483647, min3=2147483647, min4=2147483647, min5=2147483647;
        std::vector<kmer_t> minL (100,3074457345618258602);
        std::vector<kmer_t> corr_kmer (100,3074457345618258602);
    
        for(int i = 0; i<modified_sets[j].size(); i++)
        {
            for (int k = 0; k < 100; k++)
            {
                //if(j == 0 && rank == 0)
                //{   
               // printf("%d %ld %d %d %d\n",rank, modified_sets[j][i], Ax[k], Bx[k], Px[k]);
               // printf("%d %ld\n",rank, (Ax[k] * modified_sets[j][i] + Bx[k]) % Px[k]);
                //}
                /*if (((Ax[k] * modified_sets[j][i] + Bx[k]) % Px[k]) < minL[k])
                {
                    minL[k] = ((Ax[k] * modified_sets[j][i] + Bx[k]) % Px[k]);
                    corr_kmer[k] = modified_sets[j][i];
                }*/
                //oq1 = MPI_Wtime ();
                /*if ( j== 0 && i == 0 && k == 0)
                {
                    oq1 = MPI_Wtime ();
                }*/
                //val = modified_sets[j][i]*2;
                val = (((Ax[k]%Px[k])*(modified_sets[j][i]%Px[k]))%Px[k] + (Bx[k]%Px[k]))%Px[k];
                A1[j][k] = modified_sets[j][i];
                /*if ( j== 0 && i == 0 && k == 0)
                {
                    oq2 = MPI_Wtime ();
                    oq += oq2 - oq1;
                }*/
                //oq2 = MPI_Wtime ();
                //oq += oq2 - oq1;
            //if (((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]) < minL[k])
                //double tq1 = MPI_Wtime ();
                
                
                if (val <minL[k])
                {
                //minL[k] = ((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
                    //minL[k] = val;
                    //A1[j][k] = modified_sets[j][i];
                    corr_kmer[k] = modified_sets[j][i];
                }
                
                
                //double tq2 = MPI_Wtime ();
                //tq += tq2 - tq1;
                
            }
        
        }
        for (int k = 0; k < 100; k++)
        {
            A1[j][k] = corr_kmer[k];
            //if(j == 0 && rank == 0)
            //{
       // printf("%d %ld\n", rank, A1[j][k]);
            //}
    //A1[j][1] = min2;
    //A1[j][2] = min3;
    //A1[j][3] = min4;
    //A1[j][4] = min5;
        }
        corr_kmer.clear();
        corr_kmer.shrink_to_fit();
        minL.clear();
        minL.shrink_to_fit();
    }
    //printf ("%d Average time for hash across all procs (secs): %f \n", rank, oq);
    //printf ("%d Average time for update across all procs (secs): %f \n", rank, tq);
    
}

void get_hash_value_queires(std::vector<std::vector<kmer_t>> &modified_sets, kmer_t **A1)
{
    //double oq, tq, oq1, oq2 = 0.0;
    kmer_t val;
    kmer_t min_val, corr_k;
    //std::cout<<oq<< " "<<tq<<"\n";
    for(int j = 0; j<modified_sets.size(); j++)
    {
    //int min1=2147483647, min2=2147483647, min3=2147483647, min4=2147483647, min5=2147483647;
        //std::vector<kmer_t> minL (100,1074457345618108602);
        //std::vector<kmer_t> corr_kmer (100,1074457345618108602);
        for (int k = 0; k < 100; k++)
        {
            min_val = Px[k];
            for(int i = 0; i<modified_sets[j].size(); i++)
            {
                val = ((((Ax[k]%Px[k])*(modified_sets[j][i]%Px[k]))%Px[k] + (Bx[k]%Px[k]))%Px[k]);
                if (val < min_val)
                {
                    //A1[j][k] = modified_sets[j][i];
                    corr_k = modified_sets[j][i];
                    min_val = val;
                    
                }
            }
            A1[j][k] = corr_k;
        }
        /*for (int k = 0; k < 100; k++)
        {
            A1[j][k] = corr_kmer[k];
            //if(j == 0 && rank == 0)
            //{
       // printf("%d %ld\n", rank, A1[j][k]);
            //}
    //A1[j][1] = min2;
    //A1[j][2] = min3;
    //A1[j][3] = min4;
    //A1[j][4] = min5;
        }*/
        //corr_kmer.clear();
        //corr_kmer.shrink_to_fit();
        //minL.clear();
        //minL.shrink_to_fit();
    }
    //printf ("%d Average time for hash across all procs (secs): %f \n", rank, oq);
    //printf ("%d Average time for update across all procs (secs): %f \n", rank, tq);
    
}

int convert_to_int (char given_char, int len)
{
     if (given_char == 'A')
         return 0;
     else if (given_char == 'C')
         return 1;
     else if (given_char == 'G')
         return 2;
     else if (given_char == 'T')
         return 3;
     else {        
            printf (" Error: char: %d %s unrecognizable \n", len, given_char);
            exit(0);
     }
}

char convert_to_char (char given_char, int len)
{
     if (given_char == 'A')
         return 'T';
     else if (given_char == 'C')
         return 'G';
     else if (given_char == 'G')
         return 'C';
     else if (given_char == 'T')
         return 'A';
     else if (given_char == 'N')
         return 'N';
     else if (given_char == 'Y')
         return 'Y';
     else if (given_char == 'S')
         return 'S';
     else if (given_char == 'R')
         return 'R';
     else if (given_char == 'I')
         return 'I';
     else if (given_char == 'E')
         return 'E';
     else {        
            printf ("%d Error: char!: %d %s unrecognizable \n",rank, len, given_char);
            exit(0);
     }
}

void Sliding_window_queires (char *ptr, size_t length, int *num_queries,
                     std::vector<std::unordered_map<kmer_t, std::vector<int>> > Tl, int start, int total_subjects)
{
    size_t p=0;
    int max_set_length = 0;
    std::vector<kmer_t> kmer_set;
    std::vector<kmer_t> rev_set;
    std::vector<kmer_t> hash_kmers;
    std::vector<int> rev_set_tracker;
    int total_queries = 0;
    const std::string s0("/home/trahman/Asymm/OutPut/Mp1_all_mn_out_");
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
    
    //TopHit th[10] = {0};
    //std::cout<< th[2].sub<<" "<<th[2].score<<"\n";
    TopHit *th = (TopHit *)calloc(total_subjects+1, sizeof(TopHit));
    //int x;
    //std::cout<< th[node_threashold].sub<<" "<<th[node_threashold].score<< " "<< th[21].sub<<" "<<th[21].score<<"\n";
    //std::vector<int> min_val (100);
    //std::vector<kmer_t> corr_k (100);
    
    
    for(; ptr[p]!='>' && p<length; p++) { }

 
    kmer_t kmer = 0; 
    kmer_t min_kmer = 0;
    kmer_t rev_min_kmer = 0;
    kmer_t rev_kmer = 0; 
    kmer_t min_lmer_freq = 0; 
    int min_pos = 0; 
    int m_pos = 0;
    int min_val;
    int val;
    int p1 = 73;
    int p2 = 31;
    kmer_t corr_k;
    while(p<length) {
        //std::cout<<ptr[p];
        assert(ptr[p]=='>'); 
        //std::cout<<ptr[p];
    
        for(; p<length && ptr[p]!='\n'; p++) {}//std::cout<<ptr[p];}//Read the name 
        //std::cout<<"\n";
        p++; 
        total_queries++;

        if(p+100 > length) break; 

        kmer=0;
        rev_kmer=0;
        int i;
        std::string s;
        std::string str;
        int contig_len = 0;
        int tt = 0; //tracking N
        /*if (rank == 9)
        {
            printf("[%d %d]\n", rank, total_queries);
        }*/
        for(i=0; !isspace(ptr[p]) && i<100-1; i++) {
            //N Cheking
            if(tt > 0)
            {
                tt--;
            }
            if (ptr[p] == 'N' || ptr[p] == 'Y' || ptr[p] == 'S' || ptr[p] == 'R' || ptr[p] == 'I' || ptr[p] == 'E')
            {
                tt = 100;
            }
            //N Checking
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            /*
            if(rank == 0)
            {
                std::cout<<ptr[p]<<"\n";
            }*/
            s.push_back(convert_to_char(ptr[p], contig_len));
            str.push_back(ptr[p]);
            p++;
            contig_len++;
        /*    if (rank == 9 && total_queries == 1091)
            {
                if (rank == 9)
                {
                    printf("%d (---) %d\n", rank, contig_len);
                }
            } */
        }
        //int chkr = 0;
        /*
        if (rank == 9)
        {
        printf("%d {---} %d\n", rank, contig_len);
        }*/
        while(p<length && !isspace(ptr[p])) {
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //N Cheking
            if(tt > 0)
            {
                tt--;
            }
            if (ptr[p] == 'N' || ptr[p] == 'Y' || ptr[p] == 'S' || ptr[p] == 'R' || ptr[p] == 'I' || ptr[p] == 'E')
            {
                tt = 100;
            }
            //N Checking
            
            //kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            //recalculate_min_lmer(kmer, &min_kmer, &min_lmer_freq, &min_pos);
            /*if(rank == 0 && chkr == 0)
            {
                //printf("rank %d , %ld\n", rank, kmer);
                chkr += 1;
            }*/
            s.push_back(convert_to_char(ptr[p], contig_len));
            str.push_back(ptr[p]);
            //s.push_back(ptr[p]);
            p++;
            contig_len++;
            recalculate_min_kmer(str.substr(contig_len - 100, 100), &min_kmer, &min_pos, &m_pos);
            rev_set.push_back(min_kmer);
            if (min_pos > 0)
            {
                rev_set_tracker.push_back(1);
            }
            else
            {
                rev_set_tracker.push_back(0);
            }
        } 
        //std::cout << rank << " "<<s <<"\n";
        reverse(s.begin(), s.end());
        //std::cout << rank << " "<<s <<"\n";
        int tracker = 0;
        int itr = 0;
        for(i=0; i<100-1; i++) {
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //rev_kmer = kmer_shift(rev_kmer, char_to_el(s[tracker]));
            tracker++;
            //s.push_back(ptr[p]);
            //std::cout<<rank<<" "<<ptr[p]<<"\n";
            //s.push_back(convert_to_char(ptr[p]));
            //p++;
            //contig_len++;
            //s.push_back(convert_to_char(ptr[p]));
            
        }
        while(tracker<contig_len) {
            //rev_kmer = kmer_shift(rev_kmer, char_to_el(s[tracker]));
            tracker++;
            recalculate_min_kmer(s.substr(tracker - 100, 100), &rev_min_kmer, &min_pos, &m_pos);
            if(rev_min_kmer <= rev_set[contig_len-100-itr])
            {
                if(rev_set_tracker[contig_len-100-itr] == 1)
                {
                    kmer_set.push_back(rev_min_kmer);
                }
            }
            else
            {
                if(rev_set_tracker[contig_len-100-itr] == 1)
                {
                    kmer_set.push_back(rev_set[contig_len-100-itr]);
                }
            }
            
            itr++;
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            //std::cout<<rank<<" "<<ptr[p]<<"\n";
            //s.push_back(convert_to_char(ptr[p]));
            //s.push_back(ptr[p]);
            //p++;
            //contig_len++;
           // rev_set.push_back(kmer);
            //s.push_back(convert_to_char(ptr[p]));
            //s.push_back(ptr[p]);
        }
        rev_set.clear();
        rev_set.shrink_to_fit(); 
        rev_set_tracker.clear();
        rev_set_tracker.shrink_to_fit(); 
        /*
        if (rank == 0 &&  total_queries == 1)
            
            {
                std::cout <<total_queries << " " << kmer_set.size()<<"\n";
            }*/ 
        if(kmer_set.size() > 0)
        {
        /*sort(kmer_set.begin(), kmer_set.end());
        kmer_t prev=kmer_set[0];
        for(int i = 1; i < (int)(kmer_set.size()); i++)
        {
			     if (kmer_set[i] != prev) {
			         set_of_distinct_kmers.push_back(prev);
               prev=kmer_set[i];
			     }
		     
        }
                  
        set_of_distinct_kmers.push_back(prev);*/
        std::unordered_map<int, int> umap;
        int sub = 0;
        int top_hit = 0;
        for (int k = 0; k < node_threashold; k++)
        {
            //min_val = Px[k];
            min_val = Px[k];
            
            for(int i = 0; i<kmer_set.size(); i++)
            {
                //val = int(((((Ax[k]%Px[k])*(kmer_set[i]%Px[k]))%Px[k] + (Bx[k]%Px[k]))%Px[k]));
                val = int((Ax[k] * kmer_set[i] + Bx[k])%Px[k]);
                //val = (Ax[k]%Px[k]) * val;
                //std::cout<<val<<"\n";
                //val = val % Px[k];
                //int b = Bx[k] % Px[k];
                //val = (val % Px[k] + Bx[k] % Px[k]) % Px[k]; 
                //kmer_t h1 = int(((((160%p1)*(kmer_set[i]%p1))%p1 + (400%p1))%p1));
                //kmer_t h2 = int(((((1667%p2)*(kmer_set[i]%p2))%p2 + (5700%p2))%p2));
                //val = h1 + k*h2;
                if (val < min_val)
                {
                    //A1[j][k] = modified_sets[j][i];
                    corr_k = kmer_set[i];
                    min_val = val;
                    
                }
            }
            //printf("%d %d\n", rank, k);
            auto it = Tl[k].find(corr_k);
            //printf("%d %d\n", rank, k);   
            /*if (rank == 0 &&  total_queries == 1)
            
            {
                std::cout <<total_queries << " " << top_hit << " "<< sub << " "<<"\n";
            } */       
            if(it != Tl[k].end())
            {
                
                for (int elements = 0; elements < it->second.size(); elements++)
                {
                    auto iter = umap.find(it->second[elements]);
                    if (iter != umap.end())
                    {
                        iter->second = iter->second + 1;
                        /*if (top_hit < iter->second)
                        {
                            top_hit = iter->second;
                            sub = it->second[elements];
                        }*/
                    }
                    else
                    {
                        umap.insert({{it->second[elements], 1}});
                        /*
                        if (top_hit < 1)
                        {
                            top_hit = 1;
                            sub = it->second[elements];
                        }*/
                    }
                    
                    
                    
                    /* Activate this again */
                    /*
                    if (th[it->second[elements]].sub == total_queries)
                    {
                        th[it->second[elements]].score += 1;
                        if (top_hit < th[it->second[elements]].score)
                        {
                            top_hit = th[it->second[elements]].score;
                            sub = it->second[elements];
                        }
                    }
                    else
                    {
                        th[it->second[elements]].sub = total_queries;
                        th[it->second[elements]].score = 1;
                        if (top_hit < th[it->second[elements]].score)
                        {
                            top_hit = th[it->second[elements]].score;
                            sub = it->second[elements];
                        }
                        
                    }*/
                }
            }
            //hash_kmers.push_back(corr_k);
        }
        /*
        if (rank == 0 &&  total_queries == 1)
            
            {
                std::cout <<total_queries << " " << top_hit << " "<< sub << " "<<"\n";
            }*/
            
        for(auto & x:umap)
            {
                /*if (max_t < x.second)
                {
                    max_t = x.second;
                    max_sub = x.first -1;
                    unmapped = 1;
                }*/
                
                if (x.second >=3)
                {
                  fprintf(f,"%d %d %d\n",  start + total_queries - 1, x.first-1, x.second);
                }
            }
        
        
        /* Activate later
        fprintf(f,"%d %d %d\n",  start + total_queries - 1, sub - 1, top_hit);
        */
        umap.clear();
        }
        /*
        else
        {
            
            
            if (rank == 0)
            {
                std::cout<< start + total_queries - 1<<"\n";
            }
            
        }*/
        
        /*
        for(int m = 0; m<kmer_set.size(); m++)
        {    
            for(int g = 0; g < 100; g++)
            {
                min_val.push_back(Px[g]);
                corr_k.push_back(0);
            } 
            int h1 = int(((((160%p1)*(kmer_set[m]%p1))%p1 + (400%p1))%p1));
            int h2 = int(((((1667%p2)*(kmer_set[m]%p2))%p2 + (5700%p2))%p2));    
            
            for (int k =0; k<100; k++)
            {
                val = h1 + k * h2;
                if (val < min_val[k])
                {
                    //A1[j][k] = modified_sets[j][i];
                    corr_k[k] = kmer_set[m];
                    min_val[k] = val;
                    
                }
            }       
            initial_sets.push_back(corr_k);
            min_val.clear();
            min_val.shrink_to_fit();
            corr_k.clear();
            corr_k.shrink_to_fit();
        }
        } */  
               
        kmer_set.clear();
        kmer_set.shrink_to_fit();
  /*  if (max_set_length < set_of_distinct_kmers.size())
    {
        max_set_length = set_of_distinct_kmers.size();
    } */
        //initial_sets.push_back(hash_kmers);
        //hash_kmers.clear();
        //hash_kmers.shrink_to_fit();
        p++; 
        p++;
        //std::cout<<"After "<<rank<<" "<<ptr[p]<<"\n";
    }
    free(th);
  //printf("%d %d\n", rank, max_set_length);
  //int avg = 0;
  //*M_for_individual_process = max_set_length;
  *num_queries= total_queries;
  fclose(f);
  //printf("%d %d %d\n", rank, num_kmers, avg);
  //printf("total queries %d\n", total_queries);
  //num_kmers = max_set_length;
}

void Sliding_window (char *ptr, size_t length, int *M_for_individual_process, int *num_subjects,
                     std::vector<MinHashPairs> &initial_sets, int s_index)
{

    size_t p=0;
    int max_set_length = 0;
    std::vector<kmer_t> rev_set;
    std::vector<int> pos_set;
    std::vector<kmer_t> kmer_set;
    std::vector<kmer_t> kmer_set_pos;
    std::vector<kmer_t> set_of_distinct_kmers;
    std::vector<int> set_of_distinct_pos;
    std::vector<kmer_t> set_of_distinct_kmers_rev;
    std::vector<int> set_of_distinct_pos_rev;
    std::vector<kmer_t> set_of_dist_kmers;
    std::vector<MinHashPairs> set_of_dist_minhash_pairs;
    int total_subjects = 0;
  
    for(; ptr[p]!='>' && p<length; p++) { }

 
    kmer_t kmer = 0; 
    kmer_t min_kmer = 0;
    kmer_t rev_min_kmer = 0;
    kmer_t rev_kmer = 0; 
    kmer_t min_lmer_freq = 0; 
    int min_pos = 0;
    int ex_pos = 0;
    while(p<length) {
        assert(ptr[p]=='>'); 

    
        for(; p<length && ptr[p]!='\n'; p++) {
        /*if(rank == 0)
            {
                std::cout<<ptr[p];
            }*/
            }//Read the name 
            //std::cout<<"\n";
        p++; 
        total_subjects++;

        if(p+100 > length) break; 

        kmer=0;
        rev_kmer=0;
        int i;
        std::string s;
        std::string str;
        int contig_len = 0;

        for(i=0; !isspace(ptr[p]) && i<100-1; i++) {
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            //std::cout<<rank<<" "<<ptr[p]<<"\n";
            //s.push_back(ptr[p]);
            //std::cout<<rank<<" "<<ptr[p]<<"\n";
            /*if(rank == 0)
            {
                std::cout<<ptr[p]<<"\n";
            }*/
            /*if (rank == 14)
            {
              std::cout<<rank<<" "<<ptr[p]<<"\n";
            }*/
            s.push_back(convert_to_char(ptr[p], contig_len));
            str.push_back(ptr[p]);
            p++;
            contig_len++;
            //s.push_back(convert_to_char(ptr[p]));
            
        }
      //  printf("%d --- %d\n", rank, contig_len);
        while(p<length && !isspace(ptr[p])) {  
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            
            //recalculate_min_lmer(kmer, &min_kmer, &min_lmer_freq, &min_pos);
            /*if (rank == 14)
            {
              std::cout<<rank<<" "<<ptr[p]<<"\n";
            }*/
            s.push_back(convert_to_char(ptr[p], contig_len));
            str.push_back(ptr[p]);
            //s.push_back(ptr[p]);
            p++;
            contig_len++;
            //rev_set.push_back(min_kmer);
            recalculate_min_kmer(str.substr(contig_len - 100, 100), &min_kmer, &min_pos, &ex_pos);
            rev_set.push_back(min_kmer);
            pos_set.push_back(ex_pos);
            
            //s.push_back(convert_to_char(ptr[p]));
            //s.push_back(ptr[p]);
        } 
        /*if(isspace(ptr[p]))
        {
            printf("required condition\n");
        }*/
        /*
        for(int u= 0; u<rev_set.size(); u++)
        {
            printf("%ld ", rev_set[u]);
        }
        printf("End \n");*/
        //std::cout << rank << " "<<s <<"\n";
        reverse(s.begin(), s.end());
        //std::cout << rank << " "<<s <<"\n";
        int tracker = 0;
        int itr = 0;
        for(i=0; i<100-1; i++) {
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //rev_kmer = kmer_shift(rev_kmer, char_to_el(s[tracker]));
            tracker++;
            //s.push_back(ptr[p]);
            //std::cout<<rank<<" "<<ptr[p]<<"\n";
            //s.push_back(convert_to_char(ptr[p]));
            //p++;
            //contig_len++;
            //s.push_back(convert_to_char(ptr[p]));
            
        }

        while(tracker<contig_len) {
            //rev_kmer = kmer_shift(rev_kmer, char_to_el(s[tracker]));
            tracker++;
            //std::cout << tracker - 10 << "--   --"<< tracker<<"\n";
            recalculate_min_kmer(s.substr(tracker - 100, 100), &rev_min_kmer, &min_pos, &ex_pos);
            /*if (rev_min_kmer == 182726455 || rev_set[contig_len-10-itr] == 182726455)
            {
                std::cout<<tracker<<";;"<< rev_min_kmer << " " << rev_set[contig_len-10-itr] <<" " <<s.substr(tracker - 10, 10)<<"\n";
            }*/
            if(rev_min_kmer <= rev_set[contig_len-100-itr])
            {
                //if(rev_set_tracker[contig_len-10-itr] == 1)
                //{
                    if (rev_min_kmer != 3074457345618258602)
                    {
                        kmer_set.push_back(rev_min_kmer);
                        kmer_set_pos.push_back(contig_len-100-itr + (100 - 1) - (ex_pos + KMER_LENGTH-1));
                    }
               /*     if (tracker == 10)
            {
                printf("+-%ld-+\n", kmer_set[0]); 
            }*/
                //}
            }
            else
            {
                //if(rev_set_tracker[contig_len-10-itr] == 1)
                //{
                    kmer_set.push_back(rev_set[contig_len-100-itr]);
                    kmer_set_pos.push_back(pos_set[contig_len-100-itr] + contig_len-100-itr);
                //}
            }
            /*
            if (tracker == 10)
            {
                printf("--%ld--\n", kmer_set[0]); 
            }*/
            itr++;
      //kmer = (kmer_t)((kmer<<2) | (kmer_t)(convert_to_int(ptr[p]))) & (kmer_t)KMER_MASK;
            //kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            //std::cout<<rank<<" "<<ptr[p]<<"\n";
            //s.push_back(convert_to_char(ptr[p]));
            //s.push_back(ptr[p]);
            //p++;
            //contig_len++;
           // rev_set.push_back(kmer);
            //s.push_back(convert_to_char(ptr[p]));
            //s.push_back(ptr[p]);
        }
        rev_set.clear();
        rev_set.shrink_to_fit(); 
        pos_set.clear();
        pos_set.shrink_to_fit();
        /*for(int u= 0; u<kmer_set.size(); u++)
        {
           // printf("%ld ", kmer_set[u]);
        }*/
        //printf("%d %d\n", rank, set_of_distinct_kmers.size()); 
        
       // printf("End \n");
       if(contig_len >= 500)
        {
        //sort(kmer_set.begin(), kmer_set.end());
        kmer_t prev=kmer_set[0];
        //set_of_distinct_pos.push_back(kmer_set_pos[0]);
        int i ;
        for(i = 1; i < (int)(kmer_set.size()); i++)
        {
			     if (kmer_set[i] != prev) {
			         set_of_distinct_kmers.push_back(prev);
               prev=kmer_set[i];
               set_of_distinct_pos.push_back(kmer_set_pos[i-1]);
               
			     }
		     
        }
                  
        set_of_distinct_kmers.push_back(prev);
        set_of_distinct_pos.push_back(kmer_set_pos[i-1]);
        //printf("Wsize =%d\n", set_of_distinct_kmers.size());
        }
        
        for (int yu = set_of_distinct_pos.size() - 1; yu >= 0; yu --)
        {
            //printf("%d %d", set_of_distinct_pos[yu], yu);
            set_of_distinct_pos_rev.push_back(set_of_distinct_pos[yu]);
            set_of_distinct_kmers_rev.push_back(set_of_distinct_kmers[yu]);
        }
        /*
        if ( rank == 0)
        {
            for (int yu = 0; yu < set_of_distinct_pos_rev.size (); yu ++)
            {
                printf(" yu %d pos %d ", yu, set_of_distinct_pos_rev[yu]);
            }
            printf("end\n");
        }
        printf("\n Wsize =%d %d\n", set_of_distinct_kmers.size(), set_of_distinct_kmers_rev.size()); */        
        kmer_set.clear();
        kmer_set.shrink_to_fit();
        kmer_set_pos.clear();
        kmer_set_pos.shrink_to_fit();
        if (max_set_length < set_of_distinct_kmers.size())
        {
            max_set_length = set_of_distinct_kmers.size();
        }
        int j = 0;
        int k = 0;
        std::vector<kmer_t> minL (node_threashold,3074457345618258602);
        std::vector<kmer_t> corr_kmer (node_threashold,3074457345618258602);
     /*   if (rank == 0)
        {
                    printf(" k %d j %d ", k, j);
        }*/
        std::vector<std::unordered_map<kmer_t, int> > Tl(node_threashold);
        
        while (j < set_of_distinct_kmers_rev.size() && k < set_of_distinct_kmers_rev.size())
        
        {
            if ((set_of_distinct_pos_rev[k] - set_of_distinct_pos_rev[j]) <= 1000)
            {
                for(int g = 0; g < node_threashold; g++)
                {
           
                    kmer_t val = (Ax[g]*set_of_distinct_kmers_rev[k] + Bx[g]) % Px[g];
            
                    if (val <minL[g])
                    {
                //minL[k] = ((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
                        minL[g] = val;
                        corr_kmer[g] = set_of_distinct_kmers_rev[k];
                    }
                }  
                k += 1;
     
            }
            else
            {
                for(int g = 0; g < node_threashold; g++)
                {
                    //initial_sets.push_back(MinHashPairs{g, corr_kmer[g], total_subjects+s_index});
                    auto it = Tl[g].find(corr_kmer[g]);
                        
                    if(it == Tl[g].end())
                    {
                    
                        Tl[g].insert({corr_kmer[g], 1});
                    }
                    minL[g] = 3074457345618258602;
                }
                j += 1;
                k = j;
                //w_len += 1
            }
        }
        for(int g = 0; g < node_threashold; g++)
        {
            auto it = Tl[g].find(corr_kmer[g]);
                        
            if(it == Tl[g].end())
            {
                    
                Tl[g].insert({corr_kmer[g], 1});
            }
            //initial_sets.push_back(MinHashPairs{g, corr_kmer[g], total_subjects+s_index});
            //minL[g] = 1074457345618108602;
        }
        for(int g = 0; g < node_threashold; g++)
        {
        
            for(auto & x:Tl[g])
            {
                initial_sets.push_back(MinHashPairs{g, x.first, total_subjects+s_index});
            }
            Tl[g].clear();
        }
        Tl.clear();
        Tl.shrink_to_fit();
        /*
        while (j < set_of_distinct_kmers_rev.size() && k < set_of_distinct_kmers_rev.size())
        
        {
            if ((set_of_distinct_pos_rev[k] - set_of_distinct_pos_rev[j]) <= 1000)
            {
                for(int g = 0; g < 30; g++)
                {
           
                    kmer_t val = (Ax[g]*set_of_distinct_kmers_rev[k] + Bx[g]) % Px[g];
            
                    if (val <minL[g])
                    {
                //minL[k] = ((Ax[k] * kmer_sets[j][i] + Bx[k]) % Px[k]);
                        minL[g] = val;
                        corr_kmer[g] = set_of_distinct_kmers_rev[k];
                    }
                }  
                k += 1;
     
            }
            else
            {
                for(int g = 0; g < 30; g++)
                {
                    initial_sets.push_back(MinHashPairs{g, corr_kmer[g], total_subjects+s_index});
                    minL[g] = 3074457345618258602;
                }
                j += 1;
                k = j;
                //w_len += 1
            }
        }
        for(int g = 0; g < 30; g++)
        {
            initial_sets.push_back(MinHashPairs{g, corr_kmer[g], total_subjects+s_index});
            //minL[g] = 1074457345618108602;
        }*/
        /*
        for(int h = 0; h < set_of_distinct_kmers.size(); h++)
        {
            printf("%ld,", set_of_distinct_kmers[h]);
        }*/
    //    printf("\n MinHash Pair Size %d\n", set_of_dist_minhash_pairs.size());
        //initial_sets.push_back(set_of_dist_kmers);
        
        set_of_distinct_kmers.clear();
        set_of_distinct_kmers.shrink_to_fit();
        
        set_of_distinct_pos.clear();
        set_of_distinct_pos.shrink_to_fit();
        set_of_distinct_kmers_rev.clear();
        set_of_distinct_kmers_rev.shrink_to_fit();
        
        set_of_distinct_pos_rev.clear();
        set_of_distinct_pos_rev.shrink_to_fit();
        set_of_dist_kmers.clear();
        set_of_dist_kmers.shrink_to_fit();
        p++; 
        p++;
        //printf("After %c", ptr[p]);
        //std::cout<<"After "<<rank<<" "<<ptr[p]<<"\n";
    }
  //printf("Max %d %d\n", rank, max_set_length);
  //int avg = 0;
  //printf("\n MinHash Pair Size %d\n", set_of_dist_minhash_pairs.size());
  *M_for_individual_process = max_set_length;
  *num_subjects= total_subjects;
  //printf("%d %d %d\n", rank, num_kmers, avg);
  //printf("total subjects %d\n", total_subjects);
  //num_kmers = max_set_length;
}
void print_kmer_count_timers()
{

    MPI_Reduce(&sl_time, &global_sl_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for sl across all procs (secs): %f \n", 
                            (double)global_sl_time);

    MPI_Reduce(&mod, &global_mod_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for mod across all procs (secs): %f \n", 
                            (double)global_mod_time);

    MPI_Reduce(&rd, &global_rd_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for read table across all procs (secs): %f \n", 
                            (double)global_rd_time);
    
    MPI_Reduce(&hash, &global_hash_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for hash across all procs (secs): %f \n", 
                            (double)global_hash_time);
    MPI_Reduce(&fl_time, &global_fl_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for flatten across all procs (secs): %f \n", 
                            (double)global_fl_time);

    MPI_Reduce(&comm_time, &global_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for comm across all procs (secs): %f \n", 
                            (double)global_comm_time);
                            
    MPI_Reduce(&ag_time, &global_ag_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for ag across all procs (secs): %f \n", 
                            (double)global_ag_time);

    MPI_Reduce(&qsl, &global_qsl_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for qsl across all procs (secs): %f \n", 
                            (double)global_qsl_time);

    MPI_Reduce(&qhash, &global_qhash_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for qhash across all procs (secs): %f \n", 
                            (double)global_qhash_time);
    MPI_Reduce(&ss, &global_ss_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for ss across all procs (secs): %f \n", 
                            (double)global_ss_time);
    MPI_Reduce(&out, &global_out_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for out across all procs (secs): %f \n", 
                            (double)global_out_time);



}

/*
void print_kmer_count_timers()
{

    MPI_Reduce(&sl_time, &global_sl_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for sl across all procs (secs): %f \n", 
                            (double)global_sl_time/(double)size);

    MPI_Reduce(&mod, &global_mod_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for mod across all procs (secs): %f \n", 
                            (double)global_mod_time/(double)size);

    MPI_Reduce(&rd, &global_rd_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for read table across all procs (secs): %f \n", 
                            (double)global_rd_time/(double)size);
    
    MPI_Reduce(&hash, &global_hash_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for hash across all procs (secs): %f \n", 
                            (double)global_hash_time/(double)size);
    MPI_Reduce(&fl_time, &global_fl_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for flatten across all procs (secs): %f \n", 
                            (double)global_fl_time/(double)size);

    MPI_Reduce(&comm_time, &global_comm_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for comm across all procs (secs): %f \n", 
                            (double)global_comm_time/(double)size);
                            
    MPI_Reduce(&ag_time, &global_ag_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for ag across all procs (secs): %f \n", 
                            (double)global_ag_time/(double)size);

    MPI_Reduce(&qsl, &global_qsl_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for qsl across all procs (secs): %f \n", 
                            (double)global_qsl_time/(double)size);

    MPI_Reduce(&qhash, &global_qhash_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for qhash across all procs (secs): %f \n", 
                            (double)global_qhash_time/(double)size);
    MPI_Reduce(&ss, &global_ss_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for ss across all procs (secs): %f \n", 
                            (double)global_ss_time/(double)size);
    MPI_Reduce(&out, &global_out_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for out across all procs (secs): %f \n", 
                            (double)global_out_time/(double)size);



}
*/

void generate_set_of_subjects (char *read_data, size_t length, int s_index, char *r_data, size_t r_length, int start_index, int *M_final, int *num_subjects)
{
    //MPI_Barrier(MPI_COMM_WORLD);
    //double time_l5 = MPI_Wtime ();
    double time_l1 = MPI_Wtime ();
    read_array();
    //std::cout<<rank<<"\n";
    //std::vector< std::vector<kmer_t> > kmer_set_of_subjects;
    std::vector< MinHashPairs > kmer_set_of_subject;
    
    
    

    int M_for_individual_processes = 0;
    int M;
    int n_subjects;
    int total_subjects;
    
    Sliding_window (read_data, length, &M_for_individual_processes, &n_subjects, kmer_set_of_subject, s_index);
    free(read_data);
    int Array_size = kmer_set_of_subject.size();
    double time_l2 = MPI_Wtime ();
    sl_time = time_l2 - time_l1;
    
    
    //if (rank == 0) printf ("M and total subjects determined form all procceses %d %d\n", M, total_subjects);
    //*M_final = M;
    *num_subjects = total_subjects;
    
    
    std::vector<int> counts_number (size,0);
    std::vector<int> disp_array (size,0);
    double time_l3 = MPI_Wtime ();
    
    MPI_Datatype rowtype;
     MPI_Type_contiguous(sizeof(MinHashPairs), MPI_BYTE, &rowtype);
     MPI_Type_commit(&rowtype);
    MPI_Allreduce(&n_subjects, &total_subjects, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Allgather(&Array_size, 1, MPI_INT, 
            counts_number.data(), 1, MPI_INT, MPI_COMM_WORLD);
    double time_l4 = MPI_Wtime ();
    comm_time += time_l4 - time_l3;
    
    double time_l14 = MPI_Wtime ();
    //disp_node.push_ack(0);
    int total_counter = 0;
    for (int p = 0; p < size; p++)
    {
        //if (rank == 1)
        //{
   //         std::cout<<counts_number[p]<<" ";
            total_counter += counts_number[p];
        //}
    }
    //std::cout << "\n"; 
    
    std::vector< MinHashPairs > mhash_set(total_counter);
    for(int y = 1; y< size; y++)
    {
     // printf("%d ** %d\n", rank, counts_node[y]);
      disp_array[y] = disp_array[y-1] + counts_number[y-1];
    }
    double time_l24 = MPI_Wtime ();
    qsl += time_l24 - time_l14; 
    
    double time_l33 = MPI_Wtime ();
    
    MPI_Barrier(MPI_COMM_WORLD);
    int result1 = MPI_Allgatherv(kmer_set_of_subject.data(), Array_size, rowtype, 
            mhash_set.data(), counts_number.data(), disp_array.data(), rowtype, MPI_COMM_WORLD);
    //double time_l4 = MPI_Wtime ();
    //comm_time = time_l4 - time_l3;

    if (result1 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for string failed with return value: %d\n", rank, result1);
        MPI_Finalize();
        exit(2);
    }
    double time_l44 = MPI_Wtime ();
    comm_time += time_l44 - time_l33;
    
    double tq1 = MPI_Wtime ();
    kmer_set_of_subject.clear();
    kmer_set_of_subject.shrink_to_fit();
    std::vector<std::unordered_map<kmer_t, std::vector<int>> > Tl(node_threashold);
    
    for (int ita = 0; ita < mhash_set.size(); ita++)
    {
        auto it = Tl[mhash_set[ita].trail].find(mhash_set[ita].seq);    
        if(it != Tl[mhash_set[ita].trail].end())
        {
            it->second.push_back(mhash_set[ita].subject_id);
        }  
        else
        {
            std::vector<int> create_new;
            create_new.push_back(mhash_set[ita].subject_id);
            kmer_t n = mhash_set[ita].seq;
            Tl[mhash_set[ita].trail].insert({n, create_new});
            create_new.clear();
            create_new.shrink_to_fit();
                            
        }
    }
    
    mhash_set.clear();
    mhash_set.shrink_to_fit();
    double tq2 = MPI_Wtime ();
    hash = tq2 - tq1;
    
    double time_l11 = MPI_Wtime ();
    int n_queries;
    Sliding_window_queires (r_data, r_length, &n_queries, Tl, start_index, total_subjects);
    free(r_data);
    double time_l22 = MPI_Wtime ();
    qsl += time_l22 - time_l11;  
    //double time_l6 = MPI_Wtime ();
    //ss = time_l6 - time_l5; 
    print_kmer_count_timers();
    
    
    
}    

void genereate_hash_table(int M, int total_subjects, kmer_t **Ag_Hash_Table)
{
    //MPI_Allreduce(&u_ct, &g_ct, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //if (rank == 0) printf ("Average number of updates (secs): %d \n", g_ct);
    double time_l01 = MPI_Wtime ();
    
    double time_l02 = MPI_Wtime ();
    rd = time_l02 - time_l01;
    kmer_t** Prefix_table = new kmer_t*[100];
    double time_l1 = MPI_Wtime ();
    //int total_hash_functions = 100;
    for (int g = 0; g < 100; g++) {
        Prefix_table[g] = new kmer_t[M];
    }
    
    prefix_val(Prefix_table, M);
    
    int total_number_of_subs_in_p = kmer_sets.size();
    
    
    
    kmer_t** Hash_table = new kmer_t*[total_number_of_subs_in_p];
    
    int total_hash_functions = 100;
    for (int i = 0; i < total_number_of_subs_in_p; i++) {
        Hash_table[i] = new kmer_t[total_hash_functions];
    }
    
    // dynamically allocate memory of size `N` for each row
    
    //for (int ctr = 0; ctr < kmer_set_of_subjects.size(); ctr++)
    //{
    
    get_hash_value(Hash_table, M, Prefix_table);
    double time_l2 = MPI_Wtime ();
    hash = time_l2 - time_l1;
    //kmer_sets.clear();
    //kmer_sets.shrink_to_fit();
    //}
    //printf("%d %ld %ld %ld %ld %ld\n", rank, Hash_table[0][0],Hash_table[0][1], Hash_table[0][2], Hash_table[0][3], Hash_table[0][4]);
    //MPI_Barrier(MPI_COMM_WORLD);
    kmer_t* flatten_to_1d_array = (kmer_t*)malloc((total_number_of_subs_in_p*total_hash_functions) * sizeof(kmer_t));
    kmer_t* Recv_1d_array = (kmer_t*)malloc((total_subjects*total_hash_functions) * sizeof(kmer_t));
    double time_l7 = MPI_Wtime ();
    for(int ik = 0; ik < total_number_of_subs_in_p; ik++)
    {
        for(int jk = 0; jk < total_hash_functions; jk++)
        {
          //  printf("%d %ld ", rank, Hash_table[ik][jk]);
            flatten_to_1d_array[ik * total_hash_functions + jk] = Hash_table[ik][jk];
        }
        //printf("\n");
    } 
    double time_l8 = MPI_Wtime ();
    fl_time = time_l8 - time_l7;
    for( int i = 0 ; i < total_number_of_subs_in_p; i++ )
    {
        delete[] Hash_table[i]; // delete array within matrix
    }
// delete actual matrix
    delete[] Hash_table;
    int Array_size = total_number_of_subs_in_p*total_hash_functions;
    std::vector<int> counts_number (size,0);
    std::vector<int> disp_array (size,0);
    MPI_Barrier(MPI_COMM_WORLD);
    double time_l3 = MPI_Wtime ();
    MPI_Allgather(&Array_size, 1, MPI_INT, 
            counts_number.data(), 1, MPI_INT, MPI_COMM_WORLD);
    //disp_node.push_ack(0);
    for(int y = 1; y< size; y++)
    {
     // printf("%d ** %d\n", rank, counts_node[y]);
      disp_array[y] = disp_array[y-1] + counts_number[y-1];
    }
    //for(int y = 0; y< size; y++)
    //{
      //printf("%d ** %d ** %d\n", rank, counts_number[y], disp_array[y]);
      //disp_node.push_back(disp_node[y-1] + counts_node[y-1]);
    //}
    
    MPI_Barrier(MPI_COMM_WORLD);
    int result1 = MPI_Allgatherv(flatten_to_1d_array, Array_size, MPI_UINT64_T, 
            Recv_1d_array, counts_number.data(), disp_array.data(), MPI_UINT64_T, MPI_COMM_WORLD);
    double time_l4 = MPI_Wtime ();
    comm_time = time_l4 - time_l3;

    if (result1 != MPI_SUCCESS) {
        printf("rank: %d, MPI_Allgatherv for string failed with return value: %d\n", rank, result1);
        MPI_Finalize();
        exit(2);
    }
    
    /*if(rank == 1)
    {
        for(int ik = 0; ik < total_subjects*total_hash_functions; ik++)
        {
          //  printf("%ld-", Recv_1d_array[ik]);
        }
        //printf("\n");
        //printf("%ld %ld %ld \n", Recv_1d_array[0], Recv_1d_array[1], Recv_1d_array[2]);
    }*/
    //read_array();
    //double time_l5 = MPI_Wtime ();
    for(int i = 0; i < total_subjects; i++)
    {
        for(int j = 0; j < total_hash_functions; j++)
        {
            Ag_Hash_Table[i][j] = Recv_1d_array[i*total_hash_functions + j];
        }
    }
    //double time_l6 = MPI_Wtime ();
    //ag_time = time_l6 - time_l5;
    delete[] Recv_1d_array;
    /*
    for(int px = 0; px<100; px++)
    {
        //printf("%d %d\n", rank, Px[px]);
    }*/
  // MPI_Barrier(MPI_COMM_WORLD);
        
}

void generate_set_of_queries (const char *read_data, size_t length, int start_index, int total_subjects, int M, kmer_t **Ag_Hash_Table)
{
    //MPI_Barrier(MPI_COMM_WORLD);

    //std::cout<<rank<< "\n";
    //std::cout<<"start "<<start_index<<"\n";
    double tq1 = MPI_Wtime ();
    std::vector< std::vector<kmer_t> > Hash_table;
    
    int n_queries;
    
    double time_l1 = MPI_Wtime ();
    //Sliding_window_queires (read_data, length, &n_queries, Hash_table);
    double time_l2 = MPI_Wtime ();
    qsl = time_l2 - time_l1;
    //printf("%d - %d\n", rank, n_queries);
    /*
    if(rank == 0)
    {
        //printf("%ld %ld %ld %ld\n", kmer_set_of_queries[0][0], kmer_set_of_queries[0][1], kmer_set_of_queries[0][2], kmer_set_of_queries[0][3]);
        for(int j = 0; j< kmer_set_of_queries[1].size(); j++)
            {
                printf("%d %ld  ",rank, kmer_set_of_queries[1][j]);
            }
            printf("\n");
    }*/
    //MPI_Allreduce(&M_for_individual_processes, &M, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    //MPI_Allreduce(&n_subjects, &total_subjects, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    //if (rank == 0) printf ("M and total subjects determined form all procceses %d %d\n", M, total_subjects);
    //*M_final = M;
    //*num_subjects = total_subjects;
    
                         
    
    //generate_modified_set_queries(M, kmer_set_of_queries);
    //printf("%d - %d\n", rank, n_queries);
    //int total_number_of_subs_in_p = kmer_sets.size();
   /* kmer_t** Hash_table = new kmer_t*[n_queries];
    
    int total_hash_functions = 100;
    for (int i = 0; i < n_queries; i++) {
        Hash_table[i] = new kmer_t[total_hash_functions];
    }
    double time_l3 = MPI_Wtime ();
    //get_hash_value_queires(kmer_set_of_queries, Hash_table);
    double time_l4 = MPI_Wtime ();
    qhash = time_l4 - time_l3;*/
    //printf("%d - %d\n", rank, n_queries);
    //int max_match;
    //int sub;
    int total_hash_functions = 150;
    int y = 0;
    /*if(rank == 0)
    {
        //for(int i= 0; i< n_queries; i++)
        //{
            for(int j = 0; j< total_hash_functions; j++)
            {
                //printf("%d %ld  ",y++, Hash_table[1][j]);
            }
          //  printf("\n");
        //}
    }*/
    const std::string s0("/home/trahman/Asymm/OutPut/Mp1_all_mn_out_");
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
    /*
    for(int i= 0; i< n_queries; i++)
        {
           int max_match = 0;
           int sub = 0;
           //for(int j = 0; j < 1; j++)
           for(int j = 0; j < total_subjects; j++)
           {
                int match = 0;
                //printf("Printing j %d %d", rank, j);
                for(int k = 0; k< total_hash_functions; k++)
                {
                    //printf("%ld ", Hash_table[i][j]);
                    
                    if(Hash_table[i][k] == Ag_Hash_Table[j][k])
                    {
                        match += 1;
                    }
                }
                if(match >= 1)
                {
                    fprintf(f,"%d %d %d\n",  i+start_index, j, match); 
                }
                
                
           }
           //printf("Rank %d read %d matched with %d with %d\n", rank, i, sub, max_match); 
        }
    
    fclose(f);*/
    /*std::vector<std::pair<int, int>> vec1 = {{12\n, 2},
                                      {12, 1},
                                      {12, 3},
                                      {16, 3}, {16,2}, {16,3}, {16,1}, {16,-1}};*/
    //std::vector<int> indices_p(vec1.size());
    //std::vector<int> indices_s(suffix_count.size());
    //std::iota(indices_p.begin(), indices_p.end(), 0);
    //std::sort(indices_p.begin(), indices_p.end(), Comp_rev1(vec1));
    /*for (int i = 0; i < indices_p.size(); i++)
        std::cout << "[" << indices_p[i] 
             << "] ";*/
    double time_l5 = MPI_Wtime ();
    std::vector<std::unordered_map<int, int> > VM(n_queries);
    int j;
    for(int u = 0; u<total_hash_functions; u++)
    {
        std::vector<std::pair<kmer_t, int>> vec1;
        for(int o = 0; o< total_subjects; o++)
        {
            vec1.push_back({Ag_Hash_Table[o][u], o});
        }
        for(int o = 0; o<n_queries ; o++)
        {
            vec1.push_back({Hash_table[o][u], total_subjects+o});
        }
        std::vector<int> indices_p(vec1.size());
    //std::vector<int> indices_s(suffix_count.size());
        std::iota(indices_p.begin(), indices_p.end(), 0);
        std::sort(indices_p.begin(), indices_p.end(), Comp_rev1(vec1));
        
        for(int ita = 0; ita <vec1.size(); ita++)
        {
            
            kmer_t k_i = vec1[indices_p[ita]].first;
            int id_i = vec1[indices_p[ita]].second;
            if(id_i >= total_subjects && ita < vec1.size() - 1) // Boundary
            {
                j = ita + 1;
                
                kmer_t k_j = vec1[indices_p[j]].first;
                int id_j = vec1[indices_p[j]].second;
                
                while(k_j == k_i && j <vec1.size())
                {
                    
                    if(id_j < total_subjects)
                    {
                        
                        auto it = VM[id_i - total_subjects].find(id_j);
                        
                        if(it != VM[id_i - total_subjects].end())
                        {
                            it->second = it->second+1;
                        }  
                        else
                        {
                            VM[id_i - total_subjects].insert({{id_j, 1}});
                            
                        }
                    }
                    j = j + 1; //Boundary
                    if (j < vec1.size())
                    {
                        k_j = vec1[indices_p[j]].first;
                        id_j = vec1[indices_p[j]].second;
                    }
                }
            
            }
        }
        indices_p.clear();
        indices_p.shrink_to_fit();
        vec1.clear();
        vec1.shrink_to_fit();
    }
    double time_l6 = MPI_Wtime ();
    ss = time_l6 - time_l5;
    double tq2 = MPI_Wtime ();
    double tq = tq2 - tq1;
    printf ("%d Average time for l across all procs (secs): %f \n", rank, tq);
    double time_l7 = MPI_Wtime ();
    for(int yk = 0; yk <n_queries; yk++)
    {
        for(auto & x:VM[yk])
        {
            fprintf(f,"%d %d %d\n",  yk+start_index, x.first, x.second);
        }
    }
    double time_l8 = MPI_Wtime ();
    out = time_l8 - time_l7;
    fclose(f); 
    //MPI_Barrier(MPI_COMM_WORLD);
    //print_kmer_count_timers();
} 
