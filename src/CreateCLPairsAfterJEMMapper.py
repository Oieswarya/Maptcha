# Maptcha: An efficient parallel workflow for hybrid genome scaffolding

# Oieswarya Bhowmik, Tazin Rahman, Ananth Kalyanaraman

#      (oieswarya.bhowmik@wsu.edu, tazin.rahman@wsu.edu, ananth@wsu.edu)

# Washington State University

#

# **************************************************************************************************

# Copyright (c) 2024. Washington State University ("WSU"). All Rights Reserved.
# Permission to use, copy, modify, and distribute this software and its documentation
# for educational, research, and not-for-profit purposes, without fee, is hereby
# granted, provided that the above copyright notice, this paragraph and the following
# two paragraphs appear in all copies, modifications, and distributions. For
# commercial licensing opportunities, please contact The Office of Commercialization,
# WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
# commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

# IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
# OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
# THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
# ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
# OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
# **************************************************************************************************



Bench_file = open("~/JEM-Mapper/JEM_L2C.log", 'r')


Line_bench = Bench_file.readlines()

LR = {}
ct = 0
cut_off = 10
'''
for i in range (0,8000):
    Dict[i] = []'''
for line in Line_bench:
    line = line.strip()
    lt = line.split(' ')
    # if int(lt[0]) not in Li2:
    if int(lt[2]) >= cut_off:
        if int(lt[0]) not in LR:
            LR[int(lt[0])] = []
            LR[int(lt[0])].append(int(lt[1]))
        # ct += 1
        else:
            if int(lt[1]) not in LR[int(lt[0])]:
                LR[int(lt[0])].append(int(lt[1]))
            # ct+= 1





file_log_write = open('~/Maptcha/CLPairs.log', 'w')


for key in LR:
  for lr in LR[key]:
    file_log_write.write(f'{key} {int((lr))}\n')
