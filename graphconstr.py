import time

# get the start time
st = time.time()

dict_lr_set_bench = {}

file_paf_open = open('~/HM7_C_LPairsFromAsym.log', 'r')
#for line in tqdm(file_paf_open):
for line in file_paf_open:
    line = line.split('\n')[0].split(' ')
    line = [int(i) for i in line]
    if line[0] not in dict_lr_set_bench:
        dict_lr_set_bench[line[0]]  = []
        dict_lr_set_bench[line[0]].append(line[1])
    else:
        if line[1] not in dict_lr_set_bench[line[0]]:
            dict_lr_set_bench[line[0]].append(line[1])

file_w = open('~/HM7_CCPairsFromAsym.log', 'w')
list_contig = list(dict_lr_set_bench.keys())
#for c1 in tqdm(list_contig):
for c1 in list_contig:
    for c2 in list_contig[list_contig.index(c1):]:
        inter = set(dict_lr_set_bench[c2]).intersection(set(dict_lr_set_bench[c1]))
        if len(inter) > 0 and c1 != c2:
            # print(c1, c2)
            #for lr in inter:
            file_w.write(f'{c1} {c2} {len(inter)}\n')
                #file_w.write(f'{c1} {c2} {len(inter)}\n')

print("Done!!!")

et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')