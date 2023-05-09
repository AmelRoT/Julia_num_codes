#-------------------Task 1 ----------------------#

# ---------------flow rate------------------#
Q_in_L_per_h = 1100
Q_m3_per_s =(Q_in_L_per_h)*(10^-3/3600)

#-----------static head pressure -----------------#
rho = 1000;
g = 9.81
h_static_in_mmW = [366 340 255 315 339 351]
P_static_in_Pa = h_static_in_mmW*rho*g*10^-3 # 1mmW = 9.81Pa, 1mm = 10^-3m

#----------------------Area----------------------------#
d_in_m = [26 21 16 19.398 22.618 26]*10^-3 # diameter in meters
Area_all = (d_in_m).^2*pi/4

#--------------------U_measured------------------------#
U_measured= Q_m3_per_s./(Area_all) # ideal values

#-----------dynamic head pressure -----------------#

P_dynamic_in_Pa = (U_measured).^2*rho/2
h_dynamic_head = (U_measured).^2/(2*g)
h_dynamic_head_in_mmW = h_dynamic_head*1000
h_total_head_in_mmW=  h_dynamic_head_in_mmW+h_static_in_mmW
#--------------------U_calculated------------------------#
U_calculated_without_coeff = sqrt(2*(P_static_in_Pa[2]-P_static_in_Pa[1])/(rho*((Area_all[2]/Area_all[1])-1)))
C_d = U_measured[2]/U_calculated_without_coeff[1]
U_calculated = C_d*sqrt(2*(P_static_in_Pa[6]-P_static_in_Pa[5])/(rho*((Area_all[6]/Area_all[5])-1)))

U_calculated_without_coeff = sqrt(2*(P_static_in_Pa[2]-P_static_in_Pa[1])/(rho*((1-(Area_all[1]/Area_all[2])^2))))
for i = 2:length(P_static_in_Pa)-1
    U_calculated_without_coeff = [U_calculated_without_coeff sqrt(2*(P_static_in_Pa[i+1]-P_static_in_Pa[i])/(rho*((1-(Area_all[i]/Area_all[i+1])^2))))]

end

U_calculated_without_coeff =[U_calculated_without_coeff sqrt(2*(P_static_in_Pa[6]-P_static_in_Pa[5])/(rho*(((Area_all[6]/Area_all[5])^2-1))))]
cd_avgV = 0
for i = 1:length(P_static_in_Pa)
    cd_avgV =[cd_avgV U_measured[i]/U_calculated_without_coeff[i]]

end

U_calculated_withCoeff = sum(cd_avgV)/6*U_calculated_without_coeff

U_calculated_withCoeff = 0.965*U_calculated_without_coeff
#---------------------------Task 2------------------------#

#----------------------Q_calculated---------------------------------------#
Q_calculated = U_calculated_withCoeff.*Area_all
Q_measured = U_measured.*Area_all

Q_calculated_L_per_h = Q_calculated*3600*1000
Q_measured_L_per_h = Q_measured*3600*1000

using Statistics

std(U_measured)
mean(U_measured)
var(U_measured)
