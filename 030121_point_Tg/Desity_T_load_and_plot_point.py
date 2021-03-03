import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.dpi'] = 300

N_chain = 20
N_mono = 30
N_bead = 840

f_cross = 0.1
#seed_array = np.array([24305792, 64749587, 122819053, 138026712, 233001038, 270765856, 336854794, 699495782, 923177293, 941983002])
#seed_array = np.array([233001038])
seed_array = np.array([24305792, 64749587, 122819053, 138026712, 270765856, 336854794, 699495782, 923177293, 941983002])


N_cross = int(N_mono*f_cross)

N_dangle = 4
T_list = np.arange(20,102,2)
T_list = 0.01*T_list

rho_all = np.empty([len(T_list),0], float)
#rho_all = np.array([], dtype=np.int64).reshape(len(T_list),0)
T_all = np.empty([len(T_list),0], float)


#T_list = list(range(0.2, 1.02, 0.02))

for N_seed in range(0,len(seed_array)):
    seed_num = seed_array[N_seed]
    
    T_avg = np.empty((0,1), float)
    L_avg = np.empty((0,1), float)
    rho_avg = np.empty((0,1), float)
    
    for t in range(0,len(T_list)-1):
        
        if (round(100*T_list[t])%10 == 0):
            Tin = "{:.1f}".format(T_list[t])
            #print(Tin)
        else:
            Tin = "{:.2f}".format(T_list[t])
            #print(Tin)
            
        #a = round(100*T_list[t])%10
        #print(round(100*T_list[t]))
        #print(a)
        #print(Tin)
        #print(round(1.0))

        f_log = f"./f_cross={f_cross}/seed_{seed_num}/lammps_log_lj_point_T={Tin}_crosslinked_ver2.lammps"
        f = np.loadtxt(f_log, dtype=(float),max_rows=101, skiprows=245)
        
        L_tmp = np.array([2*f[:,1]])
        T_tmp = np.array([f[:,2]])
        rho_tmp = float(N_bead)/(np.average(L_tmp))**3
        L_avg = np.append(L_avg, np.average(L_tmp))
        T_avg = np.append(T_avg, np.average(T_tmp))
        rho_avg = np.append(rho_avg, rho_tmp) 
        
    # Load file for T=1 (Don't know why I cannot delete the decimal of 1.0) 
    f_log = f"./f_cross={f_cross}/seed_{seed_num}/lammps_log_lj_point_T=1_crosslinked_ver2.lammps"
    f = np.loadtxt(f_log, dtype=(float),max_rows=101, skiprows=244)
       
    L_tmp = np.array([2*f[:,1]])
    T_tmp = np.array([f[:,2]])
    #T_avg = np.average(T_tmp)
    #L_avg = np.average(L_tmp)
    #print(T_avg)
    rho_tmp = float(N_bead)/(np.average(L_tmp))**3   #x ** y	 = x to the y power
    L_avg = np.append(L_avg, np.average(L_tmp))
    T_avg = np.append(T_avg, np.average(T_tmp))   
    rho_avg = np.append(rho_avg, rho_tmp)  

    # Stack       
    
    #rho_all.append(rho_avg) 
    #T_all.append(T_avg)
    #np.concatenate((rho_all, rho_avg))
    rho_a = np.array([rho_avg]).T
    T_a = np.array([T_avg]).T
    #rho_all = np.column_stack((rho_all,rho_avg))
    
    rho_all = np.concatenate((rho_all, rho_a), axis = 1)
    T_all = np.concatenate((T_all, T_a), axis = 1)
    #np.vstack((rho_all, rho_avg))

# Averaging over seeds
T_plot = np.empty((1,0), float)
rho_plot = np.empty((1,0), float)
for i in range(0, len(T_all)):
    print("s")
    T_plot = np.append(T_plot, np.average(T_all[i,:]))
    rho_plot = np.append(rho_plot, np.average(rho_all[i,:]))

fig, (ax1) = plt.subplots(1, 1, sharex = True, figsize = (8, 6))
ax1.plot(T_plot, rho_plot ,"o", color='salmon', alpha=1, linewidth=3.5) 

ax1.set_xlabel("Temperature",fontsize=16)
ax1.set_ylabel("Density",fontsize=16)
#ax1.legend(prop={"size":12})
#ax1.set_xlim([time_data_2[1], time_data_2[-1]])
#ax1.set_ylim([0.1, 1E6])






'''
        c_atom = "./seed_{}/crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)

        c_bond = "./seed_{}/bond_crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)
        
        c_angle =  "./seed_{}/angle_crosslink_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)
       
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
    '''
'''
fig, (axes)=plt.subplots(1,1)
#aa = np.std([5, 6, 8, 9], ddof=1)
std_ratio = np.empty((0,1), float)
mean_ratio = np.empty((0,1), float)
for i in range(0, len(time_total)):
    std_ratio = np.append(std_ratio, np.std(ratio_newbond_all[i], ddof=1))
    mean_ratio = np.append(mean_ratio, np.mean(ratio_newbond_all[i]))
'''
'''
for t in range(0,len(N_bond_total)):
   ratio_newbond = np.append(ratio_newbond, (N_bond_total[t]-N_bond_total[0])/(N_chain*2))

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
''' 