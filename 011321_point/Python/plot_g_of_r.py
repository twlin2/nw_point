import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.dpi'] = 300

#from matplotlib.ticker import MaxNLocator

df = pd.read_excel('g_of_r_type2.xlsx')

#df.head(5)
#plt.plot(df['Prob=0.3'])


r_003_lmp = np.array(pd.Series(df['Bin(f=0.03)']))
g_003_lmp = np.array(pd.Series(df[' g(r) f=0.03']))

r_010_lmp = np.array(pd.Series(df['Bin(f=0.1)']))
g_010_lmp = np.array(pd.Series(df[' g(r) f=0.1']))

r_050_lmp = np.array(pd.Series(df['Bin(f=0.5)']))
g_050_lmp = np.array(pd.Series(df[' g(r) f=0.5']))



#r_side_ref = np.array(pd.Series(df['r_side']))
#r_end_ref = np.array(pd.Series(df['r_end']))
#r_cross_ref = np.array(pd.Series(df['r_cross']))


#U_side_ref = np.array(pd.Series(df['Uref(side)']))
#U_end_ref = np.array(pd.Series(df['Uref(end)']))
#U_cross_ref = np.array(pd.Series(df['Uref(cross)']))

fig,axes=plt.subplots(1,1)

axes.plot(r_003_lmp, g_003_lmp,"b-o", label="$f_{cross}$=0.03")
axes.plot(r_010_lmp, g_010_lmp,"c-.s", label="$f_{cross}$=0.11")
axes.plot(r_050_lmp, g_050_lmp,"r--^", label="$f_{cross}$=1")

'''
plt.plot(r_ss_lmp, U_ss_lmp,"r-", label="U (sphere-sphere_lmp)")
plt.plot(r_sside_lmp, U_sside_lmp,"r-.*", label="U (sphere-side_lmp)")
plt.plot(r_send_lmp, U_send_lmp,"r--", label="U (sphere-end(cross_lmp)")

plt.plot(r_pp_lmp, U_pp_lmp,"c1", label="U (point-point_lmp)")
plt.plot(r_pside_lmp, U_pside_lmp,"c-.", label="U (point-side_lmp)")
plt.plot(r_pend_lmp, U_pend_lmp,"c--", label="U (point-end_lmp)")
plt.plot(r_ps_lmp, U_ps_lmp,"c:", label="U (point-sphere_lmp)")

#plt.plot(r_side_ref, U_side_ref,"bo", label="U_side_ref")
#plt.plot(r_end_ref, U_end_ref,"ro", label="U_end_ref")
#plt.plot(r_cross_ref, U_cross_ref,"co", label="U_cross_ref")
'''
axes.set_xlabel("r")
axes.set_ylabel("g(r)")
axes.set_ylim(0,4)
axes.set_xlim(0,4)
#plt.title("The Title")
plt.legend()
axes.yaxis.set_ticks(np.linspace(0,4,5))
axes.xaxis.set_ticks(np.linspace(0,4,5))
#plt.legend(bbox_to_anchor=(1, 1.15), loc='lower center)
#plt.grid()
#plt.show()

plt.savefig("g_of_r_type2_different_fcross.png",dpi=300,format="png")


