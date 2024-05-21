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

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>

double k_gen_time=0.0, global_k_gen_time=0.0;
double k_shift_time=0.0, global_k_shift_time=0.0;
double kmer_count_time=0.0, global_kmer_count_time=0.0;
double sel_time=0.0, global_sel_time=0.0;
double alltoall_time=0.0, global_alltoall_time=0.0;
double alltoallv_time=0.0, global_alltoallv_time=0.0;
double partial_contig_time=0.0, global_partial_contig_time=0.0;
double preproc_time=0.0, global_preproc_time=0.0;
double gather_time=0.0, global_gather_time=0.0;
double barrier_time=0.0, global_barrier_time=0.0;
double sort_time=0.0, global_sort_time=0.0;
double mod = 0.0, global_mod_time = 0.0;
double sl_time = 0.0, global_sl_time = 0.0;
double hash = 0.0, global_hash_time = 0.0;
double comm_time = 0.0, global_comm_time = 0.0;
double qsl = 0.0, global_qsl_time = 0.0;
double qhash = 0.0, global_qhash_time = 0.0;
double ss = 0.0, global_ss_time = 0.0;
double out = 0.0, global_out_time = 0.0;
double rd = 0.0, global_rd_time = 0.0;
double ag_time = 0.0, global_ag_time = 0.0;
double fl_time = 0.0, global_fl_time = 0.0;



