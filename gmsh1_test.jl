
#------------------Gmsh for Julia ------------------#








#-------------------------------------------------#
r1 = 54855.48
r1_nakon = 27564.21
r2 = 66604.95
r2_nakon = 1665.3
r12_nakon = r1_nakon+r2_nakon
r12_prije = r1+r2

a1 = 148928.2
a2 = 21075.66
a1_nakon = 74834.61
a2_nakon = 526.94





#----------------------Stanje 50072532---------------------------#

s_total = 51847.84
s_total_prior = 109590.11
s_adin = 80090.11*(s_total/s_total_prior)
s_gasi = 29500*(s_total/s_total_prior)

#-------------------------Gasi oba------------------#
r1_prvi_g_dollar = 13956.66
r2_drugi_g_dollar = 895.56
r_total_dollar_g = r1_prvi_g_dollar+r2_drugi_g_dollar

r1_prvi_g_km = 25948.87
r2_drugi_g_km =1665.06
r_total_km_g = r1_prvi_g_km+r2_drugi_g_km




#--------------------------Adin oba ---------------------#
r1_prvi_a_dollar = 37891.18
r2_drugi_a_dollar = 283.38
r_total_dollar_a = r1_prvi_a_dollar+r2_drugi_a_dollar

r1_prvi_a_km = 70449.05
r2_drugi_a_km =526.94
r_total_km_a = r1_prvi_a_km+r2_drugi_a_km
#-----------------------------Total ----------------------#
r_total_dollar = r_total_dollar_a+r_total_dollar_g

r_total_dollar = r_total_km_a+r_total_km_g
