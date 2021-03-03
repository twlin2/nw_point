import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

N_chain = 20
N_mono = 30
f_cross = 0.1
N_cross = int(N_mono*f_cross)
N_dangle = 4
#tstep = 20000
box_len = 9.72000
t_max = 20000 #total time steps run in lammps
t_step = 5000 # dump file per t_step
t_num = int(t_max/t_step) + 1 # total number of frames


#seed_array = np.array([373121366, 387501919, 441758561, 530699369, 707321043, 523470439, 648841741, 887879167])
seed_array = np.array([51888245, 795596720, 34630924, 411132021, 566408070, 569965869, 842873502, 627909421, 170496349, 371137238, 799859330, 832405231])
fcross_array = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.03, 0.03, 0.03, 0.5, 0.5, 0.5])
#notype3 = 0
#a = seed_array[0]

#ratio_with_type3 = np.empty((0,1), float)
ratio_no_type3 = np.empty((0,1), float)
ratio_tri_to_self_cross = np.empty((0,1), float)
ratio_tri_to_newbond = np.empty((0,1), float)

for seed in range(0,len(seed_array)):
    seed_num = seed_array[seed]
    f_cross = fcross_array[seed]
    newbond = np.empty((0,3), int)
    bondtype2 = np.empty((0,3), int)
    bondtype1 = np.empty((0,3), int)
    num_selfcross = 0
    for t in range(0,1):
        tmp_time = t_max
        c_bond = "./seed_{}/bond_after_crosslinked_N_dangle_{}_seed_{}_fcross_{}_Nchain_{}_Nmono_{}_Rmin_1.1_Prob_0.3_Cseed_777_{}.txt".format(seed_num, N_dangle, seed_num, f_cross, N_chain, N_mono, tmp_time)
       #               bond_after_crosslinked_N_dangle_4_seed_51888245_fcross_0.1_Nchain_20_Nmono_30_Rmin_1.1_Prob_0.3_Cseed_777_20000
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
            quot = newbond[i][1]//(N_mono+N_cross*N_dangle)
            rem = newbond[i][1]%(N_mono+N_cross*N_dangle)
            bchain_1atom = quot + math.ceil((rem/(N_mono+N_cross*N_dangle)))
            quot2 = newbond[i][2]//(N_mono+N_cross*N_dangle)
            rem2 = newbond[i][2]%(N_mono+N_cross*N_dangle)
            bchain_2atom = quot2 + math.ceil((rem2/(N_mono+N_cross*N_dangle)))
            if bchain_1atom == bchain_2atom:
                num_selfcross += 1
                bond_selfcross = np.append(bond_selfcross, np.array([newbond[i]]), axis = 0)
                
    ratio_selfcross = num_selfcross/len(newbond)*100
    num_tri = 0; # number of triangular loop                
    

    ratio_no_type3 = np.append(ratio_no_type3, ratio_selfcross)
        
        
    # Count number of triangular loop when there is no type 3
    for j in range(0,len(bond_selfcross)):
        cross_tmp = bond_selfcross[j][2]
        
        for k in range(0,len(bondtype2)):
            if bondtype2[k][2]==cross_tmp-(N_dangle-1):
                if abs(bond_selfcross[j][1]- bondtype2[k][1])==1:
                    num_tri += 1
                    #print(j)
                break
        
    # calculate ratio of triangular loop formation in the self-crosslink bonds
    ratio_tri = num_tri/num_selfcross*100
    ratio_tri_to_self_cross = np.append(ratio_tri_to_self_cross, ratio_tri)
    
    ratio_tri_overall = num_tri/len(newbond)*100
    ratio_tri_to_newbond = np.append(ratio_tri_to_newbond, ratio_tri_overall)
            
        
        
            
    #print("self crosslink % = ", ratio_selfcross, "\n")
    
    #c_atom = "Coord_seed_{}_random_lb_1.300_angle_eps_10_f_reactive_0.10_N_mono_{}_N_Chains_{}_L_9.190.txt".format(seed_num, N_mono, N_chain)
            #Coord_seed_373121366_random_lb_1.300_angle_eps_10_f_reactive_0.10_N_mono_30_N_Chains_20_L_9.190
            #Coord_seed_373121366_random_lb_1.300_angle_eps_10_f_reactive_0.10_N_mono_30_N_chain_20_L_9.190.txt
    #a = np.loadtxt(c_atom, dtype = (float), skiprows=27)
    
    
#data = np.column_stack((seed_array, ratio_with_type3, ratio_no_type3))
#num_s = np.array([1, 2, 3, 4, 5])
num_s = np.array(range(1,len(seed_array)+1))

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
'''
fig, (ax1, ax2) = plt.subplots(1, 2, sharex = True, figsize = (12, 4.5))

#ax1.plot(num_s, ratio_with_type3, "b-o", label="with type 3")
ax1.plot(num_s, ratio_no_type3, "k-o", label="no type 3")
ax1.set_ylabel("Ratio of self cross-link (%)")
ax1.legend()


color = 'tab:red'
ax2.plot(num_s, ratio_tri_to_newbond, "-o", color=color, label="% in newbonds")
#ax2.plot(num_s, ratio_tri_to_self_cross, "y-o", label="no type 3")
ax2.set_ylabel("Ratio of triangular in new bonds/(%)", color=color)
ax2.tick_params(axis='y', labelcolor=color)
'''
ax3 = ax2.twinx()
color = 'tab:blue'
ax3.plot(num_s, ratio_tri_to_self_cross, "-o", color=color, label="no type 3")
ax3.set_ylabel("Ratio of triangular in self cross-link/(%)", color=color)
ax3.tick_params(axis='y', labelcolor=color)
'''
#ax2.legend(bbox_to_anchor=(1, 1.15), loc='best')
#ax2.legend(loc='upper right')


#fig.show()
fig.savefig("012121_with_seed51888245_Multi dangling beads_Self crosslink + triangular ratio.png",dpi=300,format="png")
    
print(np.average(ratio_no_type3))
print(np.std(ratio_no_type3, ddof=1))
print(np.average(ratio_tri_to_newbond))
print(np.std(ratio_tri_to_newbond, ddof=1))
print(np.average(ratio_tri_to_self_cross))
print(np.std(ratio_tri_to_self_cross, ddof=1))