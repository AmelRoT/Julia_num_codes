#----------------------Experiment 1-----------------------#

deltaH =(37.2-22.1)*10^-2
m_water = 0.605-0.110;
V_water = 0.5*10^-3;
Rho_water = m_water/V_water

P = Rho_water*9.81*deltaH

P_atm = 1.01325*10^5;
P_abs = P_atm+P


#-------------Experiment 2--------------------------#

#------------------------Diameter Section----------------#
#----------------Small spehere--------#
d_small = 16.05*10^-3
r_small = d_small/2;
#----------------Large spehere--------#
d_large = 24.44*10^-3
r_large = d_large/2;
#-----------------------End of Diameter Section------------#

#------------------------Volume Section---------------#

V_large = 4/3*(r_large^3)*pi;
V_small = 4/3*(r_small^3)*pi;
V_yogurt = 10^-3 # 1L
V_cooking_oil = 10^-3; # 1L
V_shampoo = 1*10^-3; # 1L
#----------------------End of Volume Section---------------#


#------------------Mass Section--------------#
m_large = 0.02 # kg
m_small = 0.005 # kg
m_yogurt = 0.99-0.045 # 1kg - the mass of the container
m_cooking_oil = 0.950-0.04
m_shampoo = 1.085-0.065#
#-------------------End of Mass Section----------#
#---------------------------Density Section----------------#
Rho_large = m_large/V_large;
Rho_small = m_small/V_small;
Rho_yogurt = m_yogurt/V_yogurt
Rho_cooking_oil = m_cooking_oil/V_cooking_oil
Rho_shampoo = m_shampoo/V_shampoo;
#--------------------------End of Density Section---------#

#---------------------------------Cooking Oil----------------------------#
h = 0.173
t1 = [0.24 0.31 0.33 0.27 0.25]
t2 = [0.1 0.15 0.12 0.16 0.1]
v_terminal1 = h./(t1)
v_terminal2 = h./(t2)

mju1 = 2/9*(r_small^2*9.81*(Rho_small-Rho_cooking_oil))./v_terminal1
mju2 = 2/9*(r_large^2*9.81*(Rho_large-Rho_cooking_oil))./v_terminal2

D_bottle = 9*10^-2
Re_number_cooking_oil = (Rho_cooking_oil*D_bottle*v_terminal1)/mju1

#-----------------------------Shampoo--------------------------------------#
h = 0.173
t1 = [3.11 2.98 3.01 3.06 2.92]
t2 = [1.15 1.21 1.13 1.16 1.11]
v_terminal1 = h./(t1)
v_terminal2 = h./(t2)

mju1 = 2/9*(r_small^2*9.81*(Rho_small-Rho_cooking_oil))./v_terminal1
mju2 = 2/9*(r_large^2*9.81*(Rho_large-Rho_cooking_oil))./v_terminal2

Re_number_Shampoo = (Rho_shampoo*D_bottle*v_terminal1)/mju1

#-----------------------------Yogurt ------------------------------------#
h = 0.173
t1 = [1.8 1.77 1.62 1.89 2.0]
t2 = [0.61 0.53 0.75 0.83 0.65 ]
v_terminal1 = h./(t1)
v_terminal2 = h./(t2)

mju1 = 2/9*(r_small^2*9.81*(Rho_small-Rho_cooking_oil))./v_terminal1
mju2 = 2/9*(r_large^2*9.81*(Rho_large-Rho_cooking_oil))./v_terminal2

Re_number_yogurt =  (Rho_yogurt*D_bottle*v_terminal1)/mju1

# ----------------------------Plotting the graph for oil 20-60 deg C ------#
using Plots
plotly()
T1 = 293.16
T = 293.16
for i = 2 : 40

T1 = T1+1;
T = [T (T1)];

end
b = -5;
a = -2
c = 6.74
m1 = (ones(1,length(T))*a)'+ b.*(T1./vec(T))+c.*(T1./vec(T)).^2
plot(m1)

mju = 0.356*exp.(15.7*((T1./T).-1))
plot(vec(T),vec(mju),linewidth=3,xlabel= "Temperature(T)[K]",ylabel = "Viscosity(mju)[Pa*s]",title = "Viscosity of Oil vs Temperature Graph")
