Bench_file = open("$HOME/Maptcha/TestInput/JEM_L2C.log", 'r')


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





file_log_write = open('$HOME/Maptcha/CLPairs_JEM.log', 'w')


for key in LR:
  for lr in LR[key]:
    file_log_write.write(f'{key} {int((lr))}\n')
