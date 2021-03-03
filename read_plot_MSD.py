import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from random import randint
colors = []

# generate random color array   
for i in range(10):
    colors.append('#%06X' % randint(0, 0xFFFFFF))
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = '12'




# MSD for D=2, Np=1
seed_array = np.array([18093043, 45410894, 124857597, 665368540, 413552841, 
                       704129413, 909322068, 941402059, 509901424, 5646578, 483096741, 148533427])

D_p = 2.0;
N_penetrant = 1;
AR = 1;

'''
MSD_seed_18093043_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_45410894_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_124857597_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_665368540_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_413552841_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_704129413_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_909322068_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_941402059_Dp_2.0_Np_1_AR_1_arithmetic_mixing_longer_time_new.txt

seed_509901424: 1 pene (D=2) OK
     seed_5646578: 1 pene (D=2) OK
     seed_483096741: 1 pene (D=2) OK
     seed_148533427: 1 pene (D=2) OK
'''
ratio_with_type3 = np.empty((0,1), float)
ratio_no_type3 = np.empty((0,1), float)
ratio_tri_to_self_cross = np.empty((0,1), float)
ratio_tri_to_newbond = np.empty((0,1), float)



MSD_data = np.empty((0,10001), float)

for seed in range(0, len(seed_array)):
    seed_num = seed_array[seed]
    f_name = "MSD_seed_{}_Dp_{}_Np_{}_AR_{}_arithmetic_mixing_longer_time_new.txt".format(seed_num, D_p, N_penetrant, AR)
    f = np.loadtxt(f_name, dtype=(float), skiprows=1)
    f_tmp = np.array([f[:,4]])
    if seed==0:
        MSD_data = np.append(MSD_data, f_tmp)
        time_data = np.array([f[:,0]]).T
    else:
        MSD_data = np.column_stack((MSD_data, f_tmp.T))
    #y = np.array([MSD_data[:,0]])
    #ax1.plot(time_data, ,"-o", color='royalblue', label="sphere")
    
# Calculate average
MSD_avg = np.empty((0,1), float)
a = len(MSD_data)
for i in range(0, len(MSD_data)):
    avg_tmp = np.average(MSD_data[i,:])
    MSD_avg = np.append(MSD_avg, avg_tmp)
#avg_tmp = np.average(MSD_data[])
#MSD_avg - np.append(np.average)
#print(np.average(ratio_no_type3))    
    

    
fig, (ax1) = plt.subplots(1, 1, sharex = True, figsize = (8, 6))
for seed in range(0, len(seed_array)):
    # alpha: transparency
    ax1.loglog(time_data, MSD_data[:,seed] ,"-", color='royalblue', alpha=.3,)

ax1.loglog(time_data, MSD_avg ,"-", color='royalblue', alpha=1, linewidth=3.5, label="D={}, Np={}".format(D_p, N_penetrant)) 
#ax1.plot(num_s, ratio_with_type3, "b-o", label="with type 3")
#ax1.plot(num_s, ratio_no_type3, "k-o", label="no type 3")
ax1.set_xlabel("Time",fontsize=16)
ax1.set_ylabel("MSD",fontsize=16)
ax1.legend(prop={"size":16})
ax1.set_xlim([time_data[1], time_data[-1]])
ax1.set_ylim([0.1, 1E6])



# MSD for D=1, Np=8

'''
MSD_seed_124857597_Dp_1.0_Np_8_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_665368540_Dp_1.0_Np_8_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_704129413_Dp_1.0_Np_8_AR_1_arithmetic_mixing_longer_time_new.txt
MSD_seed_941402059_Dp_1.0_Np_8_AR_1_arithmetic_mixing_longer_time_new.txt

'''
seed_array_2 = np.array([124857597, 665368540, 704129413, 941402059])
D_p_2 = 1.0;
N_penetrant_2 = 8;
AR_2 = 1;
MSD_data_2 = np.empty((0,1), float)

for seed_2 in range(0, len(seed_array_2)):
    seed_num_2 = seed_array_2[seed_2]
    #print(seed_2)
    #print(seed_num_2)
    f_name_2 = "MSD_seed_{}_Dp_{}_Np_{}_AR_{}_arithmetic_mixing_longer_time_new.txt".format(seed_num_2, D_p_2, N_penetrant_2, AR_2)
    f_2 = np.loadtxt(f_name_2, dtype=(float), skiprows=1)
    f_tmp_2 = np.array([f_2[:,4]])
    if seed_2==0:
        MSD_data_2 = np.append(MSD_data_2, f_tmp_2)
        time_data_2 = np.array([f_2[:,0]]).T
    else:
        MSD_data_2 = np.column_stack((MSD_data_2, f_tmp_2.T))
    #y = np.array([MSD_data[:,0]])
    #ax1.plot(time_data, ,"-o", color='royalblue', label="sphere")
    
# Calculate average
MSD_avg_2 = np.empty((0,1), float)

for i in range(0, len(MSD_data_2)):
    avg_tmp_2 = np.average(MSD_data_2[i,:])
    MSD_avg_2 = np.append(MSD_avg_2, avg_tmp_2)
    

for seed_2 in range(0, len(seed_array_2)):
    # alpha: transparency
    ax1.loglog(time_data_2, MSD_data_2[:,seed_2] ,"-", color='salmon', alpha=.3,)

ax1.loglog(time_data_2, MSD_avg_2 ,"-", color='salmon', alpha=1, linewidth=3.5, label="D={}, Np={}".format(D_p_2, N_penetrant_2)) 

ax1.set_xlabel("Time",fontsize=16)
ax1.set_ylabel("MSD",fontsize=16)
ax1.legend(prop={"size":12})
ax1.set_xlim([time_data_2[1], time_data_2[-1]])
ax1.set_ylim([0.1, 1E6])
'''
for seed in range(0,len(seed_array)):
    seed_num = seed_array[seed]
    
    for notype3 in range(0,1+1):
        newbond = np.empty((0,3), int)
        bondtype2 = np.empty((0,3), int)
        bondtype1 = np.empty((0,3), int)
        num_selfcross = 0
        for t in range(0,1):
            tmp_time = 40000
            if notype3==0:
                #c_atom = "soft_A_10_NH_NVT_fcross_0.1_Nchain_{}_Nmono_{}_Ncross_{}_{}.txt".format(N_chain, N_mono, N_cross, tmp_time)
                c_bond = "./seed_{}/bond_Crosslink_seed_{}_Rmin_1.1_Prob_0.3_WCA_NH_NVT_fcross_0.1_Nchain_{}_Nmono_{}_Ncross_{}_{}.txt".format(seed_num, seed_num, N_chain, N_mono, N_cross, tmp_time)
                #c_angle =  "angle_soft_A_10_NH_NVT_fcross_0.1_Nchain_{}_Nmono_{}_Ncross_{}_{}.txt".format(N_chain, N_mono, N_cross, tmp_time)
            if notype3==1:
                c_bond = "./seed_{}/bond_notype3_Crosslink_seed_{}_Rmin_1.1_Prob_0.3_WCA_NH_NVT_fcross_0.1_Nchain_{}_Nmono_{}_Ncross_{}_{}.txt".format(seed_num, seed_num, N_chain, N_mono, N_cross, tmp_time)
                
            #a = np.loadtxt(c_atom, dtype=(float), skiprows=9)
            b = np.loadtxt(c_bond, dtype=(int), skiprows=9)
            #c = np.loadtxt(c_angle, dtype=(int), skiprows=9)
            
         
            for i in range(0,len(b)):
                b_tmp = np.array([b[i,:]])
                if b[i,0]==3:
                    #b_t = b_tmp.T
                    newbond = np.append(newbond, b_tmp, axis = 0)
                elif b[i,0]==2:
                    bondtype2 = np.append(bondtype2, b_tmp, axis = 0)
                else:
                    bondtype1 = np.append(bondtype1, b_tmp, axis = 0)
                
                    
                
            bond_selfcross = np.empty((0,3), int)
            for i in range(0, len(newbond)):
                #print(newbond[i])
                quot = newbond[i][1]//(N_mono+N_cross)
                rem = newbond[i][1]%(N_mono+N_cross)
                bchain_1atom = quot + math.ceil((rem/(N_mono+N_cross)))
                quot2 = newbond[i][2]//(N_mono+N_cross)
                rem2 = newbond[i][2]%(N_mono+N_cross)
                bchain_2atom = quot2 + math.ceil((rem2/(N_mono+N_cross)))
                if bchain_1atom == bchain_2atom:
                    num_selfcross += 1
                    bond_selfcross = np.append(bond_selfcross, np.array([newbond[i]]), axis = 0)
                    
        ratio_selfcross = num_selfcross/len(newbond)*100
        num_tri = 0; # number of triangular loop                
        if notype3==0:
            ratio_with_type3 = np.append(ratio_with_type3, ratio_selfcross)
        if notype3==1:
            ratio_no_type3 = np.append(ratio_no_type3, ratio_selfcross)
            
            # Count number of triangular loop when there is no type 3
            for j in range(0,len(bond_selfcross)):
                cross_tmp = bond_selfcross[j][2]
                
                for k in range(0,len(bondtype2)):
                    if bondtype2[k][2]==cross_tmp:
                        if abs(bond_selfcross[j][1]- bondtype2[k][1])==1:
                            num_tri += 1
                            #cprint(j)
                        break
            
            # calculate ratio of triangular loop formation in the self-crosslink bonds
            ratio_tri = num_tri/num_selfcross*100
            ratio_tri_to_self_cross = np.append(ratio_tri_to_self_cross, ratio_tri)
            
            ratio_tri_overall = num_tri/len(newbond)*100
            ratio_tri_to_newbond = np.append(ratio_tri_to_newbond, ratio_tri_overall)
        
num_s = np.array(range(1,len(seed_array)+1))
'''

'''
plt.plot(num_s, ratio_with_type3, "b-o", label="with type 3")
plt.plot(num_s, ratio_no_type3, "r-o", label="no type 3")

ax = plt.gca()
ax.axes.xaxis.set_visible(False)

#plt.xlabel("time step")
plt.ylabel("Ratio of self cross-link (%)")
#plt.title("The Title")
plt.legend()
#plt.show()
plt.savefig("Self crosslink ratio.png",dpi=300,format="png")

fig, (ax1, ax2) = plt.subplots(1, 2, sharex = True, figsize = (12, 4.5))

#ax1.plot(num_s, ratio_with_type3, "b-o", label="with type 3")
ax1.plot(num_s, ratio_no_type3, "k-o", label="no type 3")
ax1.set_ylabel("Ratio of self cross-link (%)")
ax1.legend()


ax2.plot(num_s, ratio_tri_to_newbond, "c-o", label="% in newbonds")
#ax2.plot(num_s, ratio_tri_to_self_cross, "y-o", label="no type 3")
ax2.set_ylabel("Ratio of triangular in new bonds(%)")
#ax2.legend(bbox_to_anchor=(1, 1.15), loc='best')
ax2.legend(loc='upper left')



color = 'tab:red'
ax2.plot(num_s, ratio_tri_to_newbond, "-o", color=color, label="% in newbonds")
#ax2.plot(num_s, ratio_tri_to_self_cross, "y-o", label="no type 3")
ax2.set_ylabel("Ratio of triangular in new bonds/(%)", color=color)
ax2.tick_params(axis='y', labelcolor=color)


ax3 = ax2.twinx()
color = 'tab:blue'
ax3.plot(num_s, ratio_tri_to_self_cross, "-o", color=color, label="no type 3")
ax3.set_ylabel("Ratio of triangular in self cross-link/(%)", color=color)
ax3.tick_params(axis='y', labelcolor=color)

'''

