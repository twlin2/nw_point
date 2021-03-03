import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.dpi'] = 300

N_chain = 20
N_mono = 30
f_cross = 0.5
#N_cross = 1
N_cross = int(N_mono*f_cross)
N_dangle = 4
#tstep = 20000
#box_len = 9.72000
t_max = 200000 #total time steps run in lammps
t_step = 2000 # dump file per t_step
t_num = int(t_max/t_step) + 1 # total number of frames
seed_num = 526591451
tmp_time = 0;

#seed_array = np.array([795596720, 34630924, 411132021, 566408070, 569965869])
seed_array = np.array([371137238])
#seed_array = np.array([170496349, 627909421, 842873502])




for N_seed in range(0,len(seed_array)):
    seed_num = seed_array[N_seed]
    N_bond_total = np.empty((0,1), int)
    ratio_newbond = np.empty((0,1), float)
    time_total = np.empty((0,1), float)
    
    for t in range(0,t_num):
        tmp_time = t*t_step
        time_total = np.append(time_total, tmp_time)
        
        
        c_bond = "./seed_{}/bond_crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)
      #bond_Crosslink_preformed_seed_526591451_Rmin_1.1_Prob_0.3_WCA_NH_NPT_fcross_0.1_Nchain_20_Nmono_21_240000
        

        
        # Store number of bonds at every frame
        f_bond = open(c_bond, "r")
        for i in range(0,4):
            line = f_bond.readline()
            if i==3:
                N_bond_tmp = np.array(float(line))
                N_bond_total = np.append(N_bond_total, N_bond_tmp)
                
                
        
        b = np.loadtxt(c_bond, dtype=(float), skiprows=9)
       
        num_bond = np.array([range(1,len(b)+1)])
        b_out = np.concatenate((num_bond.T, b), axis=1)
        b_out_int = b_out.astype(np.int)
        
        # Another methold to stack column
        #b_out = np.column_stack((num_bond,b))
        #num_bond = np.array(range(1,len(b)+1))
        
     
        
    for step_tmp in range(0,len(N_bond_total)):
        ratio_newbond = np.append(ratio_newbond, (N_bond_total[step_tmp]-N_bond_total[0])/(N_chain*N_cross))
        ratio_tmp = np.array([ratio_newbond])
        ratio_input = ratio_tmp.T
    
    if N_seed==0:
        ratio_newbond_all = ratio_input
    #ratio_newbond_all = np.concatenate((a1, a2, ...))
    else:
        ratio_newbond_all = np.concatenate((ratio_newbond_all, ratio_input), axis=1)
    
    
#aa = np.std([5, 6, 8, 9], ddof=1)
std_ratio = np.empty((0,1), float)
mean_ratio = np.empty((0,1), float)
for i in range(0, len(time_total)):
    std_ratio = np.append(std_ratio, np.std(ratio_newbond_all[i], ddof=1))
    mean_ratio = np.append(mean_ratio, np.mean(ratio_newbond_all[i]))

'''
for t in range(0,len(N_bond_total)):
   ratio_newbond = np.append(ratio_newbond, (N_bond_total[t]-N_bond_total[0])/(N_chain*2))
  ''' 
for i in range(0, 1):
    #plt.plot(time_total, ratio_newbond_all[:,i], '-o', label="seed = {}".format(seed_array[i]))
    plt.errorbar(time_total, 100*mean_ratio, yerr=std_ratio, fmt='-o')
#plt.ylim(-0.1,1.0)
plt.xlabel("Time step")
plt.ylabel("Conversion of crosslinker (%)")
#plt.legend()

#plt.show()
#plt.savefig("Conversion of crosslinker_4 dangling_network_5 runs.png",dpi=300,format="png")
