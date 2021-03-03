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

# f=0.1
#seed_492311311 OK --> Nmax = 20
#seed_792038973 OK --> Nmax = 20
#seed_992027678 OK --> Nmax = 20
#seed_693628563
#seed_604359131

# f=0.03
#seed_988660928 OK --> Nmax = 14
#seed_354248931 OK --> Nmax = 6
#seed_726574756 OK --> Nmax = 19

seed_num = 988660928
f_cross = 0.03


num_chain = 20

# Load atom and bond information
a_name = "./seed_{}/Enp_1_Epp_1/all_after_crosslinked_N_dangle_4_seed_{}_Rmin_1.1_Prob_0.3_WCA_NH_NPT_fcross_{}_Nchain_20_Nmono_30_Enp_1_Epp_1_100000.txt".format(seed_num, seed_num, f_cross)
a = np.loadtxt(a_name, dtype=(float), skiprows=9)
b_name = "./seed_{}/Enp_1_Epp_1/bond_after_crosslinked_N_dangle_4_seed_{}_Rmin_1.1_Prob_0.3_WCA_NH_NPT_fcross_{}_Nchain_20_Nmono_30_Enp_1_Epp_1_100000.txt".format(seed_num, seed_num, f_cross)
b = np.loadtxt(b_name, dtype=(int), skiprows=9)


a_ID = np.array([a[:,0]]).astype(np.int)
a_mol = np.array([a[:,1]]).astype(np.int)
a_type = np.array([a[:,2]]).astype(np.int)
a_list = np.concatenate((a_ID.T, a_mol.T, a_type.T), axis=1)

# Sort atoms by its molecule number
list = []
for i in range(0, num_chain):
    n_tmp = np.empty((0,3), int)
    for j in range(0, len(a_list)):
        if (a_list[j,1] == i+1):
            a_tmp = np.array([a_list[j,:]])
            n_tmp = np.append(n_tmp, a_tmp, axis=0)
    list.append(n_tmp)

# Store bond type 3 (new bonds formed in crosslinking)
b_type3_org = np.empty((0,3), int)
for i in range(0,len(b)):
    b_tmp = np.array([b[i,:]])
    if b[i,0]==3:
        b_type3_org = np.append(b_type3_org, b_tmp, axis = 0)
        

# Index check for bond
b_type3 = b_type3_org
b_check = np.ones((len(b_type3),3), int)





# Start of counting network size

chain_check = np.empty((0,1), int)
nw_size_list = np.empty((0,1), int)

#nw_size = 1


while (len(chain_check) < num_chain):
    nw_size = 1
    print("Start counting chain size")
    
    # Delete type3 bonds which have been checked
    b_type3 = b_type3_org*b_check

    # Find the chain to start 
    # this will need to be modified to choose type3 bonds that "haven't been visited"
    for i in range(0, len(b_type3)):
        if b_type3[i,1] != 0: 
            b_a1_tmp = b_type3[i,1]
            b_a2_tmp = b_type3[i,2]
            break
        
    for i in range(0, len(a_list)):
        if (a_list[i,0] == b_a1_tmp) or (a_list[i,0] == b_a2_tmp):
            if a_list[i,2] != 4: # choose the non-crosslinker bead to find its chain number
                c_tmp = a_list[i,1]
                break
    
    # Store number of chain that has been checked
    chain_check = np.append(chain_check, c_tmp) 
    
    # Store atom info for selected chain
    # We don't need this step as we already sort atom by their molecular number
    '''
    a_tmp_chain = np.empty((0,3), int)
    for i in range(0, len(a_list)):
        if (a_list[i,1] == c_tmp):
            a_tmp = np.array([a_list[i,:]])
            a_tmp_chain = np.append(a_tmp_chain, a_tmp, axis = 0)
    '''
    
    
    stack = []   
    
    a_tmp_chain = list[c_tmp-1]
    for i in range(0, len(b_type3)):
        b_a1 = b_type3[i,1]
        b_a2 = b_type3[i,2]
        for j in range(0,len(a_tmp_chain)):
            # Use two if to prevent storing atom info for self-crosslink bond
            if (a_tmp_chain[j,0] == b_a1):
                stack.append(a_tmp_chain[j,0])
                break
            elif (a_tmp_chain[j,0] == b_a2):
                stack.append(a_tmp_chain[j,0])
                break    
            
    while (stack != []): 
        end_index = 0
        #flag = int(0)
        #print(stack)
        
        ap = stack.pop()   #
        #print(ap)
        
        # pick the first type3 bond and find the atom next to this bond
        # --> find its molecule number
        for i in range(0, len(b_type3)):
            if (b_type3[i,1] == ap) or (b_type3[i,2] == ap):
                if b_type3[i,1] == ap:
                    a_next_p = b_type3[i,2]
                    break
                elif b_type3[i,2] == ap:
                    a_next_p = b_type3[i,1]
                    break
        
        
        for i in range(0, len(b_type3)):
            for j in range(1,3):
                if (b_type3[i,j] == ap):
                    b_check[i,j] = 0
                    #flag = 1

                if (b_type3[i,j] == a_next_p):
                    b_check[i,j] = 0
                    #flag = 2

            
        # Delete the row of bond that has been checked
        #b_type3 = np.delete(b_type3, i, axis = 0)       
        #print (b_type3)
                  
        
        # Find a_next_p mol ID
        for i in range(0, len(a_list)):
            if (a_list[i,0] == a_next_p):
                c_tmp = a_list[i,1]
                break         
            
        for i in range(0, len(chain_check)):
            if chain_check[i] == c_tmp:
                # connect back to network itself --> end this branch
                #print("Crosslink back to network itself. End this branch")
                end_index = 2
                
                break
            
        #print("end_index = {}".format(end_index))
        
        if (end_index != 2):    # If this chain hasn't been checked  
            # if this chain hasn't been visited, then find the type3 bond on this chain and store the atoms number with type3 bond
            #print("nw_size before = {}".format(nw_size))
            nw_size += 1
            #print("nw_size after = {}".format(nw_size))
            chain_check = np.append(chain_check, c_tmp) 
            
            a_tmp_chain = list[c_tmp-1]
            for j in range(0, len(b_type3)):
                b_a1 = b_type3[j,1]
                b_a2 = b_type3[j,2]
                for k in range(0,len(a_tmp_chain)):
                    # Use "two if" to prevent storing atom info for self-crosslink bond
                    if (a_tmp_chain[k,0] == b_a1) and (a_tmp_chain[k,0] != a_next_p): # exclude the atom connecting the 'previou' chain
                        stack.append(a_tmp_chain[k,0])
                        break
                    elif (a_tmp_chain[k,0] == b_a2) and (a_tmp_chain[k,0] != a_next_p):
                        stack.append(a_tmp_chain[k,0])
                        break
        else: 
            print("Crosslink back to network itself. Go back to search another brance")

# Print out the maximum network size
        # ROUND 2
        
        '''
        ap = stack.pop()   #
        
        # pick the first type3 bond and find the atom next to this bond
        # --> find its molecule number
        for i in range(0, len(b_type3)):
            if (b_type3[i,1] == ap) or (b_type3[i,2] == ap):
                if b_type3[i,1] == ap:
                    a_next_p = b_type3[i,2]
                    break
                elif b_type3[i,2] == ap:
                    a_next_p = b_type3[i,1]
                    break
        
        flag = int(0)
        for i in range(0, len(b_type3)):
            for j in range(1,3):
                if (b_type3[i,j] == ap):
                    b_check[i,j] = 0
                    flag = 1
                    break
            if (flag ==  1):
                break
            
        # Delete the row of bond that has been checked
        b_type3 = np.delete(b_type3, i, axis = 0)            
                  
        
        for i in range(0, len(a_list)):
            if (a_list[i,0] == a_next_p):
                c_tmp = a_list[i,1]
                break       
            
        for i in range(0, len(chain_check)):
            if chain_check[i] == c_tmp:
                # connect back to network itself --> end this branch
                print("Crosslink back to network itself. End this branch")
                end_index = 2
                break
            
            
        if (end_index != 2):
            # if this chain hasn't been visited, then find the type3 bond on this chain and store the atoms number with type3 bond
            print("nw_size before = {}".format(nw_size))
            nw_size += 1
            print("nw_size after = {}".format(nw_size))
            chain_check = np.append(chain_check, c_tmp) 
            
            a_tmp_chain = list[c_tmp-1]
            for j in range(0, len(b_type3)):
                b_a1 = b_type3[j,1]
                b_a2 = b_type3[j,2]
                for k in range(0,len(a_tmp_chain)):
                    # Use "two if" to prevent storing atom info for self-crosslink bond
                    if (a_tmp_chain[k,0] == b_a1) and (a_tmp_chain[k,0] != a_next_p): # exclude the atom connecting the 'previou' chain
                        stack.append(a_tmp_chain[k,0])
                        break
                    elif (a_tmp_chain[k,0] == b_a2) and (a_tmp_chain[k,0] != a_next_p):
                        stack.append(a_tmp_chain[k,0])
                        break
        '''

    nw_size_list = np.append(nw_size_list, nw_size)

# Print out maximum network size
print("Maximum Network Size = {}".format(np.max(nw_size_list)))



'''
# MSD for D=2, Np=1
seed_array = np.array([18093043, 45410894, 124857597, 665368540, 413552841, 
                       704129413, 909322068, 941402059, 509901424, 5646578, 483096741, 148533427])

D_p = 2.0;
N_penetrant = 1;
AR = 1;


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