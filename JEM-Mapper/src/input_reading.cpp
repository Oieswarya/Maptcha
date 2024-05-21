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


double file_open_time=0.0, global_file_open_time=0.0;
double file_read_time=0.0, global_file_read_time=0.0;
double file_copy_time=0.0, global_file_copy_time=0.0;
double file_close_time=0.0, global_file_close_time=0.0;
double input_process_time=0.0, global_input_process_time=0.0;


char* DivideReads(MPI_File *in, const int rank, const int size, const int overlap,
                 uint64_t *nlines, size_t *data_size) {
    MPI_Offset filesize;
    MPI_Offset localsize, ori_localsize;
    MPI_Offset start;
    MPI_Offset end;
    char *chunk;
    uint64_t pad_width=1000000000; // 1GB

#ifdef DEBUG
    char proc_id[3];
    char output_file_name[25];

    sprintf(proc_id, "%d", rank); 
    strcpy(output_file_name,"out_");
    strcpy(&output_file_name[strlen(output_file_name)],proc_id);
    strcpy(&output_file_name[strlen(output_file_name)],".log");

    char debug_file_name[25];
    strcpy(debug_file_name,"debug_p");
    strcpy(&debug_file_name[strlen(debug_file_name)],proc_id);
    strcpy(&debug_file_name[strlen(debug_file_name)],".log");

    FILE *f = fopen(output_file_name, "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    FILE *f0 = fopen(debug_file_name, "w");
    if (f0 == NULL)
    {
        printf("Error opening p0 file!\n");
        exit(1);
    }
#endif

    /* figure out who reads what */

    double t3 = MPI_Wtime ();
    MPI_File_get_size(*in, &filesize);
    localsize = filesize/size;
    start = rank * localsize;
    end   = start + localsize - 1;
 
#ifdef DEBUG_IO
    printf("Proc:%d size of file is: %lld, start: %lld, end: %lld \n", rank, filesize, start, end);
#endif
    /* add overlap to the end of everyone's chunk... */
    end += overlap;

    /* except the last processor, of course */
    if (rank == size-1) end = (filesize-1);

    localsize =  end - start + 1;

    if (rank == 0) fprintf (stderr, "filesize: %llu, localsize: %llu\n", filesize, localsize);

#ifdef DEBUG_IO
    printf("Proc:%d, size of file is: %lld, start: %lld, end: %lld, localsize: %lld \n", rank, filesize, start, end, localsize);
#endif
    /* allocate memory */
    chunk = (char*) malloc( (localsize + 1)*sizeof(char));

    //create contiguous derived data type
    MPI_Datatype rowtype;
    MPI_Type_contiguous(pad_width, MPI_BYTE, &rowtype);
    MPI_Type_commit(&rowtype);

    /* Check the size of localsize */
    ori_localsize=localsize;
    size_t dest_offset=0;
    if (localsize > pad_width)
    {
      int count_size = floor(localsize/pad_width);

      if (rank == 0) fprintf (stderr, "count_size: %d\n", count_size);

      /* everyone reads in their part */
      MPI_File_read_at_all(*in, start, chunk, count_size, rowtype, MPI_STATUS_IGNORE);
      //MPI_File_read_at_all(*in, start, chunk, localsize, MPI_BYTE, MPI_STATUS_IGNORE);
      localsize -= (count_size*pad_width);
      
      assert (ori_localsize == ((count_size*pad_width) + localsize));
      start += count_size*pad_width;
      dest_offset = count_size*pad_width;
    }
    //MPI_Offset new_offset = count_size*pad_width;
    //MPI_File_read_at_all(*in, new_offset, &chunk[new_offset], localsize, MPI_BYTE, MPI_STATUS_IGNORE);
    
    MPI_File_read_at_all(*in, start, &chunk[dest_offset], localsize, MPI_BYTE, MPI_STATUS_IGNORE);
    chunk[ori_localsize] = '\0';

    if (rank == 0) fprintf (stderr, "dest_offset: %lu, localsize: %llu\n", dest_offset, localsize);

    // free datatype
    MPI_Type_free(&rowtype);

    localsize = ori_localsize;
    if (rank == 0) fprintf (stderr, "dest_offset: %lu, localsize: %llu\n", dest_offset, localsize);
    double t4 = MPI_Wtime ();
    file_read_time = t4 - t3;
    MPI_Reduce(&file_read_time, &global_file_read_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_Read (secs): %f \n", 
                            (double)global_file_read_time/(double)size);

    /*
     * everyone calculate what their actual start and end positions on a read
     * from the first '>' after start to the first '>' after the
     * overlap region starts (eg, after end - overlap + 1)
     */

#ifdef DEBUG
    /* Print data */
    for (int i=0; i<localsize; i++)
         fprintf(f0,"%c", chunk[i]);
    printf("\n");
        
    fclose(f0);

    printf("Before: Proc:%d locstart: %lld, locend: %lld, localsize: %lld, overlap: %d \n", rank, locstart, locend, localsize, overlap);
#endif

    double t11 = MPI_Wtime();
    uint64_t locstart=0, locend=localsize;
    //unsigned long long int locstart=0, locend=localsize;

    if (rank == 0) fprintf (stderr, "initial locstart: %lu, locend: %lu\n", locstart, locend);

    if (rank != 0) {
        //std::cout<<rank << " "<<locstart<<"\n";
        while(chunk[locstart] != '>') 
        {
            locstart++;
            //std::cout<<rank << " "<<locstart<<"\n";   
        }
    }
    //std::cout<<"Check "<<rank << " "<<locstart<< " " << locend<<"\n";  
    if (rank != size-1) {
        locend-=overlap;
        //std::cout<<rank << " "<<locend<<"\n";
        while(chunk[locend] != '>')
        {  
           locend++;
          // if (rank == 0)
          // {
          // std::cout<<rank << " "<<locend<<"\n"; 
          // std::cout<<rank << " "<<locend<< " " << chunk[locend]<<"\n"; 
          // }
        }
    }
    //std::cout<<"Check "<<rank << " "<<locstart<< " " << locend<<"\n"; 
    localsize = locend-locstart;
    
    if (rank == 0) fprintf (stderr, "new localsize: %llu, locstart: %lu, locend: %lu\n", localsize, locstart, locend);

#ifdef DEBUG
    printf("After: Proc:%d locstart: %lld, locend: %lld, localsize: %lld \n", rank, locstart, locend, localsize);
#endif

    /* Now let's copy our actual data over into a new array, with no overlaps */
    char *data = (char *)malloc((localsize+1)*sizeof(char));
    memcpy(data, &(chunk[locstart]), localsize);
    data[localsize] = '\0';
    free(chunk);

    /* Now we'll count the number of lines */
    *nlines = 0;
    for (int i=0; i<localsize; i++)
        if (data[i] == '>') (*nlines)++;

    *data_size = localsize;
#ifdef DEBUG
    /* Print data to corresponding output files */
    for (int i=0; i<localsize; i++)
         fprintf(f,"%c", data[i]);
    //printf("\n");

    fclose(f);
#endif

    double t12 = MPI_Wtime();
    file_copy_time = t12-t11;
    MPI_Reduce(&file_copy_time, &global_file_copy_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_copy (secs): %f \n", 
                            (double)global_file_copy_time/(double)size);

    //free(data);

    return data;
}


input_read_data perform_input_reading (const int rank, const int size,
                                       std::string &fileName, int read_length)
{

    MPI_File in;
    int ierr;
   
    const int overlap=read_length+100;
    uint64_t nlines=0, global_nlines=0, localsum;
    size_t read_data_size=0;
    input_read_data input_rdata;

    double start_t = MPI_Wtime ();

    double t1 = MPI_Wtime ();
    ierr = MPI_File_open(MPI_COMM_WORLD, fileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
    if (ierr) {
        if (rank == 0) fprintf(stderr, "Couldn't open the input FASTA file\n");
        MPI_Finalize();
        exit(2);
    }

    double t2 = MPI_Wtime ();
    file_open_time = t2 - t1;

    MPI_Reduce(&file_open_time, &global_file_open_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_Open (secs): %f \n", 
                            (double)global_file_open_time/(double)size);

    if (rank==0)
        fprintf(stderr, "Start Reading the input dataset\n");
        
    input_rdata.read_data = DivideReads(&in, rank, size, overlap, &nlines, &read_data_size);
    input_rdata.read_data_size = read_data_size;
    MPI_Scan(&nlines, &localsum, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);  
    input_rdata.start_index =  localsum-nlines;

#ifdef DEBUG_OUT
    printf("Rank %d has %d lines\n", rank, nlines);
#endif

    double t5 = MPI_Wtime ();
    MPI_File_close(&in);

    MPI_Barrier(MPI_COMM_WORLD);
    double t6 = MPI_Wtime ();
    file_close_time = t6 - t5;
    MPI_Reduce(&file_close_time, &global_file_close_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for MPI_FIle_Close (secs): %f \n", 
                            (double)global_file_close_time/(double)size);

    double end_t = MPI_Wtime ();
    input_process_time = end_t - start_t;

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Reduce(&input_process_time, &global_input_process_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for reading and storing the Reads in each proc's memory (secs): %f \n", 
                            (double)global_input_process_time);

    MPI_Reduce (&nlines, &global_nlines, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Total number of reads: %lu\n", global_nlines);
    if (rank == 0) fprintf (stderr, "Total number of reads: %lu\n", global_nlines);

    if (rank==0)
        fprintf(stderr, "Completed Reading the input dataset\n");

    return input_rdata;

}

