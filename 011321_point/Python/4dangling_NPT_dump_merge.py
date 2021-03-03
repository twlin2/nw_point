import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.dpi'] = 300


N_chain = 20
N_mono = 30

#f_cross = 0.5
#seed_array = np.array([371137238, 799859330, 832405231])

#f_cross = 0.1
#seed_array = np.array([34630924, 51888245, 411132021, 566408070, 569965869, 795596720])

f_cross = 0.5
#seed_array = np.array([170496349, 627909421, 842873502])
#seed_array = np.array([371137238, 799859330, 832405231])
seed_array = np.array([371137238])


N_cross = int(N_mono*f_cross)
#N_cross = 1
N_dangle = 4
#tstep = 20000
#box_len = 9.72000
t_max = 200000 #total time steps run in lammps
t_step = 2000 # dump file per t_step
t_num = int(t_max/t_step) + 1 # total number of frames
#seed_num = 526591451
tmp_time = 0;



for N_seed in range(0,len(seed_array)):
    seed_num = seed_array[N_seed]
    N_bond_total = np.empty((0,1), int)
    ratio_newbond = np.empty((0,1), float)
    time_total = np.empty((0,1), float)
    
    for t in range(0,t_num):
        tmp_time = t*t_step
        time_total = np.append(time_total, tmp_time)
        
        c_atom = "./seed_{}/crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)

        #crosslink_N_dangle_4_seed_371137238_fcross_0.5_Nchain_20_Nmono_30_Rmin_1.1_Prob_0.3_Cseed_777_104000
        #Crosslink_N_dangle_4_seed_371137238_fcross_0.5_Nchain_20_Nmono_30_Rmin_1.1_Prob_0.3_Cseed_777_1000.txt'
        #c_atom = str("soft_A_10_NH_NVT_fcross_0.1_Nchain_2_Nmono_20_Ncross_2_0.txt")
        c_bond = "./seed_{}/bond_crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)
        #bond_Crosslink_preformed_seed_526591451_Rmin_1.1_Prob_0.3_WCA_NH_NPT_fcross_0.1_Nchain_20_Nmono_21_240000
        
        c_angle =  "./seed_{}/angle_crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)
        
        #angle_after_crosslinked_N_dangle_4_seed_371137238_fcross_0.5_Nchain_20_Nmono_30_Rmin_1.1_Prob_0.3_Cseed_777_9000
        #Crosslink_seed_859194152_Rmin_1.1_Prob_0.3_WCA_NH_NVT_fcross_0.1_Nchain_2_Nmono_30_Ncross_3_0
        #bond_Crosslink_seed_859194152_Rmin_1.1_Prob_0.3_WCA_NH_NVT_fcross_0.1_Nchain_2_Nmono_30_Ncross_3_0
    
        # Store new box dimensions
        box_len = np.empty((0,1), float)
        f = open(c_atom,"r")
        for i in range(0,8):
            line = f.readline()
            #print(line)
            if i>=5:
                aa = np.array(str.split(line))
                len_i = float(aa[1])
                box_len = np.append(box_len, len_i)
        f.close()
        
        # Store number of bonds at every frame
        f_bond = open(c_bond, "r")
        for i in range(0,4):
            line = f_bond.readline()
            if i==3:
                N_bond_tmp = np.array(float(line))
                N_bond_total = np.append(N_bond_total, N_bond_tmp)
                
                
        a = np.loadtxt(c_atom, dtype=(float), skiprows=9)
        b = np.loadtxt(c_bond, dtype=(float), skiprows=9)
        c = np.loadtxt(c_angle, dtype=(float), skiprows=9)
        
        a_ID = a[:,0]
        a_ID_int = a_ID.astype(np.int)
        a_molID = a[:,1]
        a_molID_int = a_molID.astype(np.int)
        a_type = a[:,2]
        a_type_int = a_type.astype(np.int)
        a_charge = a[:,3]
        a_position = a[:,4:8]
        
        num_bond = np.array([range(1,len(b)+1)])
        b_out = np.concatenate((num_bond.T, b), axis=1)
        b_out_int = b_out.astype(np.int)
        
        # Another methold to stack column
        #b_out = np.column_stack((num_bond,b))
        #num_bond = np.array(range(1,len(b)+1))
        
        num_angle = np.array([range(1,len(c)+1)])
        c_out = np.concatenate((num_angle.T, c), axis=1)
        c_out_int = c_out.astype(np.int)
        
        fil = open("./seed_{}/dump__seed_{}_{}".format(seed_num, seed_num, tmp_time,) + ".txt", "w")
        
        fil.write(str("#Polymer network - time step = ") + "{}\n\n".format(tmp_time))
        fil.write("{}".format(len(a)) + str(" atoms\n"))
        fil.write("{}".format(len(b)) + str(" bonds\n"))
        fil.write("{}".format(len(c)) + str(" angles\n"))
        fil.write("{}".format(0) + str(" dihedrals\n"))
        fil.write("{}".format(0) + str(" impropers\n\n"))
        
        fil.write("{}".format(4) + str(" atom types\n"))
        fil.write("{}".format(3) + str(" bond types\n"))
        fil.write("{}".format(3) + str(" angle types\n\n"))
        
        fil.write("{:} {}".format(-box_len[0], box_len[0]) + str(" xlo xhi\n"))
        fil.write("{} {}".format(-box_len[1], box_len[1]) + str(" ylo yhi\n"))
        fil.write("{} {}".format(-box_len[2], box_len[2]) + str(" zlo zhi\n\n"))
        
        
        fil.write(str("Masses\n\n"))
        for i in range(0,4):
            fil.write("{} {}\n".format(i+1, 1.000000))
        
        fil.write(str("\nAtoms\n\n"))
        for i in range(0, len(a)):
            fil.write(str(a_ID_int[i]) + " ")
            fil.write(str(a_molID_int[i]) + " ")
            fil.write(str(a_type_int[i]) + " ")
            fil.write(str(a_charge[i]) + " ")
            fil.write(str(a_position[i])[1:-1] + "\n")
            
        
        fil.write(str("\n\nBonds\n\n"))
        for i in range (0, len(b_out_int)):
            #np.savetxt(fil, b_out_int[i,:])
            fil.write(str(b_out_int[i])[1:-1] + '\n')
        
        fil.write(str("\n\nAngles\n\n"))
        for i in range (0, len(c_out_int)):
            #np.savetxt(fil, b_out_int[i,:])
            fil.write(str(c_out_int[i])[1:-1] + '\n')
        
        fil.write("\n")
        
        #fil1.write(str(all_salts[b]) + "\t{0:.6f}\n".format(concentrations[b]))
        
        fil.close()
        
    for step_tmp in range(0,len(N_bond_total)):
        ratio_newbond = np.append(ratio_newbond, (N_bond_total[step_tmp]-N_bond_total[0])/(N_chain*N_cross))
        ratio_tmp = np.array([ratio_newbond])
        ratio_input = ratio_tmp.T
    
    if N_seed==0:
        ratio_newbond_all = ratio_input
    #ratio_newbond_all = np.concatenate((a1, a2, ...))
    else:
        ratio_newbond_all = np.concatenate((ratio_newbond_all, ratio_input), axis=1)
    
    
fig, (axes)=plt.subplots(1,1)
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
    axes.errorbar(time_total, mean_ratio, yerr=std_ratio, fmt='-o')
axes.set_ylim(-0.1,1.1)
#axes.set_xlim(-1000,21000)
#axes.xaxis.set_ticks(np.linspace(0,20000,5))
axes.set_xlabel("Time step")
axes.set_ylabel("Conversion of crosslinker")
#plt.legend()

#plt.show()
#plt.savefig("Conversion of crosslinker_preformed network_5 runs.png",dpi=300,format="png")
