
# question 38
mju_glc = 1.5
u = 5.5
y = 6*10^-3
tau = mju_glc *(u)/y;
rho = 1260
Re = (rho)*u*y/mju_glc

#1.47
rho = 880;
v = 0.003
mju = v*rho;
rout = 6.02/2*10^-2
rin =  6.00/2*10^-2
l = 40*10^-2
Ffr = mju*((0.4)/(rout-rin))*l*pi*rin*2

#1.41
rout = 6.04/2*10^-2
rin =  6.00/2*10^-2
mju_oil = 0.86
u = (30)*(rout-rin)/(pi*(6*10^-2)*40*10^-2*mju_oil)

y1 = 0.011/2*(1/(1.97*10^-5))*2*3*10^-2/(10.8*pi)

y2 = acos(y1)

yf = y2*2*3*10^-2/pi

# chapter 2
# 2.10
rho_hg = 13500
rho_w = 1000
rho_oil =  888
Pair = 60*10^3-9.81*(rho_hg*0.2+rho_w+rho_oil*1.5)

## 2.12

h1 = 1000*0.18/(898)
h = h1-0.12

# 2.13

h = (14715-11772+15696)/7848
