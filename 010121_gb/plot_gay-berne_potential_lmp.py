import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['figure.dpi'] = 300


df = pd.read_excel('py_input_gb_potential_lmp_022121_compare.xlsx')

#df.head(5)
#plt.plot(df['Prob=0.3'])


r_side_lmp = np.array(pd.Series(df['r_side_lmp']))
r_end_lmp = np.array(pd.Series(df['r_end_lmp']))
r_cross_lmp = np.array(pd.Series(df['r_cross_lmp']))
r_T_lmp = np.array(pd.Series(df['r_T_lmp']))

#r_ss = np.array(pd.Series(df['rss']))
#r_se = np.array(pd.Series(df['rse']))
#r_spsp = np.array(pd.Series(df['rspsp']))


#r_ss_lmp = np.array(pd.Series(df['r(sphere-sphere)']))
#r_sside_lmp = np.array(pd.Series(df['r(sphere-side)']))
#r_send_lmp = np.array(pd.Series(df['r(sphere-end)']))

#r_pp_lmp = np.array(pd.Series(df['r(point-point)']))
#r_pside_lmp = np.array(pd.Series(df['r(point-side)']))
#r_pend_lmp = np.array(pd.Series(df['r(point-end)']))
#r_ps_lmp = np.array(pd.Series(df['r(point-sphere)']))


r_side_ref = np.array(pd.Series(df['r_side']))
r_end_ref = np.array(pd.Series(df['r_end']))
r_cross_ref = np.array(pd.Series(df['r_cross']))
r_T_ref = np.array(pd.Series(df['r_T']))

U_side_lmp = np.array(pd.Series(df['U(side-to-side)']))
U_end_lmp = np.array(pd.Series(df['U(end-to-end)']))
U_cross_lmp = np.array(pd.Series(df['U(cross)']))
#U_T_lmp = np.array(pd.Series(df['U(T)']))

#U_ss_lmp = np.array(pd.Series(df['U(sphere-sphere)']))
#U_sside_lmp = np.array(pd.Series(df['U(sphere-side)']))
#U_send_lmp = np.array(pd.Series(df['U(sphere-end)']))

#U_ss = np.array(pd.Series(df['Uss']))
#U_se = np.array(pd.Series(df['Use']))
#U_spsp = np.array(pd.Series(df['Uspsp']))




'''
U_pp_lmp = np.array(pd.Series(df['U(point-point)']))
U_pside_lmp = np.array(pd.Series(df['U(point-side)']))
U_pend_lmp = np.array(pd.Series(df['U(point-end)']))
U_ps_lmp = np.array(pd.Series(df['U(point-sphere)']))
'''

U_side_ref = np.array(pd.Series(df['Uref(side)']))
U_end_ref = np.array(pd.Series(df['Uref(end)']))
U_cross_ref = np.array(pd.Series(df['Uref(cross)']))
#U_T_ref = np.array(pd.Series(df['Uref(T)']))


plt.plot(r_side_lmp, U_side_lmp,"b-", label="U (side_lmp)")
plt.plot(r_end_lmp, U_end_lmp,"b-.", label="U (end_lmp)")
plt.plot(r_cross_lmp, U_cross_lmp,"b--.", label="U (cross_lmp)")
#plt.plot(r_T_lmp, U_T_lmp,"ko", label="U (T_lmp)")

#plt.plot(r_ss_lmp, U_ss_lmp,"b-", label="U (sphere-sphere_lmp)")
#plt.plot(r_sside_lmp, U_sside_lmp,"b-o", label="U (sphere-side_lmp)")
#plt.plot(r_send_lmp, U_send_lmp,"b--", label="U (sphere-end(cross_lmp)")

#plt.plot(r_ss, U_ss,"r-o", label="U (sphere-side_lmp)")
#plt.plot(r_se, U_se,"r-*", label="U (sphere-end_lmp)")
#plt.plot(r_spsp, U_spsp,"r-", label="U (sphere-sphere(cross_lmp)")






'''
plt.plot(r_pp_lmp, U_pp_lmp,"c1", label="U (point-point_lmp)")
plt.plot(r_pside_lmp, U_pside_lmp,"c-.", label="U (point-side_lmp)")
plt.plot(r_pend_lmp, U_pend_lmp,"c--", label="U (point-end_lmp)")
plt.plot(r_ps_lmp, U_ps_lmp,"c:", label="U (point-sphere_lmp)")
'''
plt.plot(r_side_ref, U_side_ref,"r-", label="U_side_ref")
plt.plot(r_end_ref, U_end_ref,"r-.", label="U_end_ref")
plt.plot(r_cross_ref, U_cross_ref,"r--", label="U_cross_ref")
#plt.plot(r_T_ref, U_T_ref,"k-", label="U_T_ref")

plt.xlabel("r(center to center)")
plt.ylabel("U")
plt.ylim(-5,0.5)
#plt.title("The Title")
#plt.legend()
#plt.legend(bbox_to_anchor=(1, 1.15), loc='lower center)
#plt.grid()
#plt.show()

#plt.savefig("gb_potential_lmp_022121_GB2521_finizte size only.png",dpi=300,format="png")


