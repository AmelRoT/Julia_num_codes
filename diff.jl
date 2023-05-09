#---------------------Differential Equations-----------------


#------------------------Euler Method -------------------------------------
using Plots
plotly()


a = 1
b = 1.5
n = 100
h = (b-a)/n

xi = a:h:b #interval
yi = 1 #intial conditions
#f(xi,yi) = (2/xi)*yi+xi.^2*exp.(xi)
f(xi,yi) = 2*xi
s = yi; #keeping track of yi+1 values
for i in 1:length(xi)-1 # n = lenght(xi)-1

yf = yi+h*f.(xi[i],yi);
yi = yf;
s = [s yf];

end


plot(xi,vec(s),linewidth=3)
scatter!(xi,vec(s))
f2(x) = exp(x)
plot!(xi, f2.(xi),linewidth = 3)
f1 = x->70/9*exp.(-0.3*x)-43/9*exp.(-1.2*x)
f1 = x->x.^2+2*x+2-2*(x+1)*log.(x+1)
f1 = x->sin.(x)+cos.(x)
f1 =x-> 1/2*exp.(3*(cos.(x)-1))

f1 = x-> (exp.(x-1))/(2*x)
f1 = x-> x.^2*(exp.(x)-exp(1))


#---------------------------Heun's Method----------------------------------


a = 0
b = 6
n = 15
h = (b-a)/n

xi = a:h:b #interval
yi = 1 #intial conditions
yn1 = yi;
k = -1
f(xi,yi) = k*yi+(1-k)*cos.(xi)-(1+k)*sin.(xi)
s = yi; #keeping track of yi+1 values
for i in 1:length(xi)-1 # n = lenght(xi)-1

yf = yi+h*f.(xi[i],yi)
yi = yf;
s = [s yf]
end
s1 = s[1];
for i in 1:n

ynew = yn1+(f.(xi[i],yn1)+f.(xi[i+1],s[i+1]))*h/2 #s starts from a value and it is corrected by i+1
yn1 = ynew
s1 = [s1 ynew]

end


scatter!(xi,vec(s1),linewidth=3)

scatter(xi,vec(s1),linewidth=3)


#--------------Midpoint Method------------------------



a = 0
b = 6
n = 15
h = (b-a)/n

xi = a:h:b #interval
yi = 1 #intial conditions
k = -1
f(xi,yi) = k*yi+(1-k)*cos.(xi)-(1+k)*sin.(xi)

s = yi; #keeping track of yi+1 values
s1 = yi
yn1 = yi;
for i in 1:n
ym = yn1+h/2*f.(xi[i],yn1)
ynew = yn1+h*(f.((xi[i]+h/2),ym)) #s starts from a value and it is corrected by i+1
yn1 = ynew
s1 = [s1 ynew]

end
scatter!(xi,vec(s1),linewidth=3)

#---------------------------Backward Euler Method-----------------
#bounds a---------b
using SymPy
using Plots
plotly()
x = symbols("x")
y = symbols("y")
a = 0
b = 1
n = 2 #number of intervals
h = (b-a)/n

xi = a:h:b #interval
yi = 1/2 #intial conditions
k = -1
i = 1;
g(xi,yi) = -yi*log.(yi)
y =ynew-> -ynew+yi+h*g.(xi[i+1],ynew)
y1 = derivative(y)

s = yi; #keeping track of yi+1 values
s1 = yi;
for i in 1:n
s = NewtonMethod.(y,y1,yi)
s1 = [s1 s]
yi = s;

end
plot(xi,exp.(xi),linewidth=3)
scatter!(xi,vec(s1),linewidth=3)

#---------------------------Trapezoid Method-----------------
#bounds a---------b
using Plots
plotly()
a = 0;
b = 1;
n = 1
h =(b-a)/n #length of the interval sections
yi = 1/2 #intial yi
xi = a:h:b #whole inerval of x
s = yi
f(xi,yi) =-yi*log.(yi)

#intial conditions
k = -1;
i = 1;
g(xi,yi) =-yi*log.(yi)
y =ynew-> -ynew+yi+h/2*(g.(xi[i],yi)+g.(xi[i+1],ynew))
y1 = derivative(y)

s = yi; #keeping track of yi+1 values
s1 = yi;
for i in 1:n
s = NewtonMethod.(y,y1,yi)
s1 = [s1 s]
yi = s;

end
plot(xi,exp.(xi),linewidth=3)
scatter!(xi,vec(s1),linewidth=3)


#-----------------------------Runge Kutta 4th order--------

a = 0;
b = 1;
n = 100
h =(b-a)/n #length of the interval sections
yi = 0 #intial yi
xi = a:h:b #whole inerval of x
s = yi
f1(xi,yi) = exp.(-xi.^2)

for i in 1:n

k1 = f1.(xi[i],yi)
k2 = f1.((xi[i]+h/2),(yi+k1/2*h))
k3 = f1.((xi[i]+h/2),(yi+k2/2*h))
k4 = f1.((xi[i]+h),(yi+k3*h))

k = 1/6*(k1+2*k2+2*k3+k4)
ynew = yi+k*h;
yi = ynew;
s =[s ynew]
end
y =x-> 70/9*exp(-0.3*x)-43/9*exp(-1.2*x);
scatter!(xi,vec(s),linewidth=3)
plot(xi,y.(xi),linewidth = 3)




#--------------------------Runge Kutta 2nd order ---------------
#---for plots -------------
using Plots
plotly()
#---------------------------------------
#Modified Euler Method - Runge Kutta Method c2 = 1/2
a = 0;
b = 1;
n = 100;
h =(b-a)/n #length of the interval sections
yi = 1 #intial yi
xi = a:h:b #whole inerval of x
s = yi
f(xi,yi) =(2-2*xi*yi)*(xi.^2+1).^-1
c2 =1/2;
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = yi;
for i in 1:n

k1 = f.(xi[i],yi);
k2 = f((xi[i]+alfa*h),(yi+beta*k1*h))

ynew = yi+h*(c1*k1+c2*k2);
yi = ynew;
s = [s ynew]
end
y =x-> sqrt(x^2+2*x+6)-1
scatter!(xi,vec(s),linewidth=3)
plot(xi,y.(xi),linewidth = 3)
#----------------------------------------
#Midpoint Runge Kutta Method c2 = 1
a = 1;
b = 2;
n = 10;
h = (b-a)/n;
yi = 0;
xi = a:h:b;
f(xi,yi) = xi+3*yi/xi
c2 = 1;
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = yi;
for i in 1:n

k1 = f.(xi[i],yi);
k2 = f((xi[i]+alfa*h),(yi+beta*k1*h))

ynew = yi+h*(c1*k1+c2*k2);
yi = ynew;
s = [s ynew]
end

#------------------------------------------
#Heun's Runge Kutta Method c2 = 3/4
a = 1;
b = 2;
n = 10;
h = (b-a)/n;
yi = 0;
xi = a:h:b;
f(xi,yi) = xi+3*yi/xi
c2 = 3/4;
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = yi;
for i in 1:n

k1 = f.(xi[i],yi);
k2 = f((xi[i]+alfa*h),(yi+beta*k1*h))

ynew = yi+h*(c1*k1+c2*k2);
yi = ynew;
s = [s ynew]
end
#-------------------------------------------------

#y1 =x-> 25/16*exp.(-4/3*x)+3/4*x-9/16
#plot(xi,y1.(xi),linewidth=3)

c2 = 3/4;
c2 = 1;

f(xi,yi) = -xi*yi+4*xi/yi
y = x-> (4-3*exp(-x.^2)).^(1/2)
#------------------------

y2 = 1;
for i in 1:length(xi)-1

y1g = sum((-1).^i*(factorial(big(i))).^-1*xi.^(2*i))
y2 = [y2 y1g]

end

y4 =x-> 1-x.^2+(x.^4)/2-(x.^6)/6+(x.^6)/(factorial(4))-(x.^10)/(factorial(5))+(x.^12)/(factorial(6))-(x.^14)/(factorial(7)+(x.^16)/(factorial(8)))-(x.^18)/(factorial(9))



plot(xi,y4.(xi),linewidth=3)

#--------------------------Runge Kutta 3rd order -calssical ---------------

a = 0;
b = 1.5;
n = 10;
h = (b-a)/n;
yi = 3;
xi = a:h:b;
s = yi
f1(xi,yi) =-1.2*yi+7*exp.(-0.3*xi)

for i in 1:n

k1 = f1.(xi[i],yi)
k2 = f1.((xi[i]+h/2),(yi+k1/2*h))
k3 = f1.((xi[i]+h),(yi-k1*h+2*k2*h))


k = 1/6*(k1+4*k2+k3)
ynew = yi+k*h;
yi = ynew;
s =[s ynew]

end
scatter!(xi,vec(s),linewidth=3)


#---------------Higher order, systems  - Runge Kutta method -----------------------

using Plots;
plotly()



using SymPy;

a = 0;
b = 2;
n = 100
h =(b-a)/n #length of the interval sections
xi = a:h:b # whole inerval of x

z1 = 0;
z2 = 0;
zi = [z1;z2];
A = [0 1;-3/2 -1/2]
f1(z) = A*z+[0;5]
s = zi;
for i in 1:n
f1.(z1,z2);
zf = zi+h*f1.(z1,z2)
zi =zf;
z1 = zi[1];
z2= zi[2];
s = [s zi];
end
s1 = s[1];
s2 = s[2]
for i in 3:2:length(s)

s1 = [s1 s[i]]
end
for i in 4:2:length(s)

s2 = [s2 s[i]]
end


y =x-> 20/sqrt(23)*exp.(-x/4).*sin.(sqrt(23)/4*x)
plot(xi,y.(xi),linewidth=3)
scatter!(xi,vec(s1),linewidth=3)
scatter!(xi,vec(s2),linewidth=3)



a = 0;
b = 2;
n = 100;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x

z1 = 0;
z2 = 0;
zi = [z1;z2];
f(z1,z2) = [0 1;-3 -3]*zi+[0;1]
s = zi;
c2 = 1;
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = zi;
for i in 1:n

k1 = f.(z1,z2);
k2 = f((xi[i]+alfa*h),(zi+beta*k1*h))
zf = zi.+h*(c1*k1.+c2*k2);
zi = zf;
s = [s zf]
end

s1 = s[1];
s2 = s[2]
for i in 3:2:length(s)

s1 = [s1 s[i]]
end
for i in 4:2:length(s)

s2 = [s2 s[i]]
end


y =x-> 2/sqrt(3)*sin.(sqrt(3)/2*x)*exp.(-3/2*x)
y = x->1/3-1/3*exp.(-3/2*x)*cos.(sqrt(3)/2*x)-exp.(-3/2*x)*1/sqrt(3)*sin(sqrt(3)/2*x)
plot(xi,y.(xi),linewidth=3)

scatter!(xi,vec(s1),linewidth=3)
scatter!(xi,vec(s2),linewidth=3)





a = 0;
b = 100;
n = 300;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
b = 0.5;
z1 = pi/2;
z2 = 0;
zi = [z1;z2];
f(z1,z2) = [z2 -9.81/3*sin.(z1)-b/(3)*z2]
s = zi;
for i in 1:n

zf = zi+h*vec(f.(z1,z2))
zi =zf;
z1 = zi[1];
z2= zi[2];
s = [s zi];
end
s1 = s[1];
s2 = s[2]
for i in 3:2:length(s)

s1 = [s1 s[i]]
end
for i in 4:2:length(s)

s2 = [s2 s[i]]
end


y =x-> 2/sqrt(3)*sin.(sqrt(3)/2*x)*exp.(-3/2*x)
y = x->1/3-1/3*exp.(-3/2*x)*cos.(sqrt(3)/2*x)-exp.(-3/2*x)*1/sqrt(3)*sin(sqrt(3)/2*x)
plot(xi,y.(xi),linewidth=3)

scatter!(xi,vec(s1),linewidth=3)
scatter(xi,vec(s2),linewidth=3)


#------------RK 2nd order works ---------------------
a = 0;
b = 100;
n = 5000;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
b = 0.5;
z1 = pi/2;
z2 = 0;
zi = [z1;z2];
f(z1,z2) = [z2 -9.81/3*sin.(z1)]
f(z1,z2) = [z2 -9.81/3*sin.(z1)-b/(3)*z2]

c2 = 1;
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = zi;
for i in 1:n

k1 = f.(z1,z2);
k2 = f.((z1+beta*k1[1]*h),(z2+beta*k1[2]*h))
zf = zi+vec(h*(c1*k1.+c2*k2));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];

end
s1 = s[1];
s2 = s[2]
for i in 3:2:length(s)

s1 = [s1 s[i]]
end
for i in 4:2:length(s)

s2 = [s2 s[i]]
end
s3 = s[1,:]
s4 = s[2,:]


s1 = s[1,:]
s2 = s[2,:]

plot(xi,vec(s1),linewidth=3)
p2 = plot(xi,vec(s1),linewidth=3)
p1 = plot(xi,vec(s2),linewidth=3)
p3 = plot(xi,vec(s3),linewidth=3)
p4 = plot(xi,vec(s4),linewidth=3)

plot(p1,p4,layout=(2,1))
plot(p2,p3,layout=(2,1))


#--------------Works - tested with the book ---------------------------------
a = 0;
b = 1;
n = 20
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 = 1
z2 = 1
zi = [z1;z2]
s = zi;
c2 = 1/2
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
f(xi,z1,z2) = [3 2;4 1]*[z1;z2]+[-(2*xi.^2+1)*exp.(2*xi);((xi.^2)+2*xi-4)*exp.(2*xi)]
for i in 1:n

k1 = f.(xi[i],z1,z2);
k2 = f.((xi[i]+alfa*h),(z1+beta*k1[1]*h),(z2+beta*k1[2]*h))
zf = zi+vec((h*(c1*k1+c2*k2)));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];

end

s1 = s[1,:];
s2 = s[2,:];


#same here
a = 1;
b = 2;
n = 100
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 = 2
z2 = 8
z3 = 6
zi = [z1;z2;z3]
c2 = 1/2
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = zi
f(xi,z1,z2,z3) = [z2;z3;(-1/xi)*z3+2/((xi).^2)*z2-2/((xi).^3)*z1+8-2/(xi).^3]
for i in 1:n

k1 = f.(xi[i],z1,z2,z3);
k2 = f.((xi[i]+alfa*h),(z1+beta*k1[1]*h),(z2+beta*k1[2]*h),(z3+beta*k1[3]*h))
zf = zi+vec(h*(c1*k1.+c2*k2));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];
z3 = zi[3];

end

s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];



#last problem

a = 0;
b = 3;
n = 15
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 = 3
z2 = -1
z3 = 9
zi = [z1;z2;z3]
c2 = 1/2
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = zi
f(xi,z1,z2,z3) = [z2;z3;z3/xi-(3/xi.^2)*z2+4/(xi.^3)*z1+5*log.(xi)+9]
for i in 1:n

k1 = f.(xi[i],z1,z2,z3);
k2 = f.((xi[i]+alfa*h),(z1+beta*k1[1]*h),(z2+beta*k1[2]*h),(z3+beta*k1[3]*h))
zf = zi+vec(h*(c1*k1.+c2*k2));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];
z3 = zi[3];

end
s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];












y = x->2*x-x.^(-1)+x.^2+x.^3-1

plot(xi,y.(xi),linewidth=3)
scatter!(xi,vec(s1),linewidth=3)




#same here x2
a = 0;
b = 1;
n = 10
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 = 3
z2 = -1
z3 = 1
zi = [z1;z2;z3]
c2 = 1/2
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = zi
f1(xi,z1,z2) = [1 2 -2;0 1 1;1 2 0]*[z1;z2;z3]+[exp.(-xi);-2*exp(-xi);exp.(-xi)]
for i in 1:n

k1 = f.(xi[i],z1,z2,z3);
k2 = f.((xi[i]+alfa*h),(z1+beta*k1[1]*h),(z2+beta*k1[2]*h),(z3+beta*k1[3]*h))
zf = zi+vec(h*(c1*k1.+c2*k2));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];
z3 = zi[3];
end

s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];


y1 = x-> -3*exp.(-x)-3*sin.(x)+6*cos.(x)


y = x->2*x-x.^(-1)+x.^2+x.^3-1
y1 = x->2*exp.(3*x)+3*exp.(-2*x)+x
y2 = x->-8*exp.(-2*x)+exp.(4*x)-2*exp.(3*x)+sin.(x)
y3 = x->2*exp.(4*x)-4*exp.(3*x)-exp(-2*x)-2

plot(xi,y1.(xi),linewidth=3)
plot(xi,y2.(xi),linewidth=3)
plot(xi,y3.(xi),linewidth=3)

scatter!(xi,vec(s1),linewidth=3)
scatter!(xi,vec(s2),linewidth=3)
scatter!(xi,vec(s3),linewidth=3)













#4order RK
a = 0;
b = 2;
n = 15
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 = 3
z2 = -1
z3 = 9
zi = [z1;z2;z3]
s = zi;
f1(xi,z1,z2,z3) = [0 1 0;0 0 1;4 4 -1]*[z1;z2;z3]+[0;0;0]

for i in 1:n

k1 = f1.(xi[i],z1,z2,z3)
k2 = f1.((xi[i]+h/2),(z1+k1[1]/2*h),(z2+k1[2]/2*h),(z3+k1[3]/2*h))
k3 = f1.((xi[i]+h/2),(z1+k2[1]/2*h),(z2+k2[2]/2*h),(z3+k2[3]/2*h))
k4 = f1.((xi[i]+h),(z1+k3[1]*h),(z2+k3[2]*h),(z3+k3[3]*h))
k = 1/6*(k1+2*k2+2*k3+k4)
zf = zi+k*h;
zi = zf;
s =[s zf]
z1 = zi[1];
z2 = zi[2];
z3 = zi[3];
end

s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];


y1 = xi->exp(-xi)+exp.(2*xi)+exp(-2*xi)
y2 = xi->-8*exp(-2*xi)+exp.(4*xi)-2*exp.(3*xi)+sin.(xi)
y3 = xi->2*exp(4*xi)-4*exp(3*xi)-3*exp(-2*xi)-2

plot(xi,y1.(xi),linewidth=3)
plot!(xi,y2.(xi),linewidth=3)
plot!(xi,y3.(xi),linewidth=3)
scatter!(xi,s1)
scatter!(xi,s2)
scatter!(xi,s3)


using Plots
plotly()

# this one works 4th order RK method - with the book
a = 0;
b = 25;
n = 1000
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 =500
z2 = 3000
bg = 1.1
bl = 0.00025
dg = 0.0005
dl = 0.7
zi = [z1;z2]
s = zi;
f1(xi,z1,z2) = [bl*z1*z2-dl*z1;bg*z2-dg*z1*z2]

for i in 1:n

k1 = f1.(xi[i],z1,z2)
k2 = f1.((xi[i]+h/2),(z1+k1[1]/2*h),(z2+k1[2]/2*h))
k3 = f1.((xi[i]+h/2),(z1+k2[1]/2*h),(z2+k2[2]/2*h))
k4 = f1.((xi[i]+h),(z1+k3[1]*h),(z2+k3[2]*h))
k = 1/6*(k1+2*k2+2*k3+k4)
zf = zi+k*h;
zi = zf;
s =[s zf]
z1 = zi[1];
z2 = zi[2];
end

s1 = s[1,:];
s2 = s[2,:];

scatter(xi,vec(s1))
scatter!(xi,vec(s2))

y = x->-1/2*exp.(2*x)+x.^2+2*x-1/2
y = x->1/2*exp.(2*x)+x.^2-1/2

#new one for the 4RK from with 3 eq. ## ---exam problem ----
a = 1;
b = 2;
n = 10;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x
z1 = 2
z2 = 8
z3 = 6
zi = [z1;z2;z3]
s = zi;
f1(xi,z1,z2,z3) = [0 1 0;0 0 1;-2/(xi.^3) 2/(xi.^2) -1/(xi)]*[z1;z2;z3]+[0;0;8-2/(xi.^3)]

for i in 1:n

k1 = f1.(xi[i],z1,z2,z3)
k2 = f1.((xi[i]+h/2),(z1+k1[1]/2*h),(z2+k1[2]/2*h),(z3+k1[3]/2*h))
k3 = f1.((xi[i]+h/2),(z1+k2[1]/2*h),(z2+k2[2]/2*h),(z3+k2[3]/2*h))
k4 = f1.((xi[i]+h),(z1+k3[1]*h),(z2+k3[2]*h),(z3+k3[3]*h))
k = 1/6*(k1+2*k2+2*k3+k4)
zf = zi+k*h;
zi = zf;
s =[s zf]
z1 = zi[1];
z2 = zi[2];
z3 = zi[3];
end

s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];


y = x->2*exp.(3*x)+3*exp.(-2*x)+x



using Plots
plotly();
a = 0;
b = 4.5;
n = 2000
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x

z1 = 0.355028053887817;
z2 = -0.2588194503792807;
zi = [z1;z2]

s = zi;
f1(xi,z1,z2) = [z2;xi*z1]

for i in 1:n

k1 = f1.(xi[i],z1,z2)
k2 = f1.((xi[i]+h/2),(z1+k1[1]/2*h),(z2+k1[2]/2*h))
k3 = f1.((xi[i]+h/2),(z1+k2[1]/2*h),(z2+k2[2]/2*h))
k4 = f1.((xi[i]+h),(z1+k3[1]*h),(z2+k3[2]*h))
k = 1/6*(k1+2*k2+2*k3+k4)
zf = zi+k*h;
zi = zf;
s =[s zf]
z1 = zi[1];
z2 = zi[2];
end

s1 = s[1,:];
s2 = s[2,:];

scatter!(xi,s1)
scatter!(xi,s2)




a = 0;
b = 2*pi;
n = 300;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x

z1 = 1
z2 = 0
z3 = 0
z4 = 1
zi = [z1;z2;z3;z4]
s = zi;
f1(xi,z1,z2,z3,z4) = [z3;z4;-z1*(z1.^2+z2.^2).^(-3/2);-z2*(z1.^2+z2.^2).^(-3/2)]

for i in 1:n

k1 = f1.(xi[i],z1,z2,z3,z4)
k2 = f1.((xi[i]+h/2),(z1+k1[1]/2*h),(z2+k1[2]/2*h),(z3+k1[3]/2*h),(z4+k1[4]/2*h))
k3 = f1.((xi[i]+h/2),(z1+k2[1]/2*h),(z2+k2[2]/2*h),(z3+k2[3]/2*h),(z4+k2[4]/2*h))
k4 = f1.((xi[i]+h),(z1+k3[1]*h),(z2+k3[2]*h),(z3+k3[3]*h),(z4+k3[4]*h))
k = 1/6*(k1+2*k2+2*k3+k4)
zf = zi+k*h;
zi = zf;
s =[s zf]
z1 = zi[1];
z2 = zi[2];
z3 = zi[3];
z4 = zi[4];
end

s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];
s4 = s[4,:];

using Plots
plotly()


scatter(xi,s1)
scatter!(xi,s2)

scatter(xi,s3)
scatter!(xi,s4)


y = x->exp.(-2*x)-6*x*exp.(-2*x);
y1 = x->-exp.(-2*x)+2*x*exp.(-2*x);

plot(xi,y.(xi),linewidth=3);
plot!(xi,y1.(xi),linewidth=3)


n1 =x-> cos.(x)+sin.(x)
n2 = x-> 2*cos.(x)

plot(xi,n1.(xi),linewidth=4)
plot!(xi,n2.(xi),linewidth=4)

#-------------
a = 0;
b = 100;
n = 6000;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x

z1 = 1
z2 = 1
z3 = 1
zi = [z1;z2;z3]
s = zi
c2 = 1;
c1 = 1-c2;
alfa = 1/(2*c2)
beta = alfa;
s = zi;
f(xi,z1,z2,z3) = [0 1 0;0 0 1;-1 -3 -3]*zi+[0;0;-4*sin(xi)]
for i in 1:n

k1 = f.(xi[i],z1,z2,z3);
k2 = f.(xi[i],(z1+beta*k1[1]*h),(z2+beta*k1[2]*h),(z3+beta*k1[3]*h))
zf = zi+vec(h*(c1*k1.+c2*k2));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];
z3 = zi[3];
end


s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];

scatter(xi,s1)
scatter!(xi,s2)
scatter!(xi,s3)




#----------------------------------------------
a = 0;
b = 2;
n = 100;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x

z1 = 3
z2 = 0.2
zi = [z1;z2]
s = zi
f(xi,z1,z2) = [(-z1+z2)*exp.(1-xi)+0.5*z1;z1-z2.^2]
for i in 1:n

k1 = f.(xi[i],z1,z2);
k2 = f.(xi[i],(z1+beta*k1[1]*h),(z2+beta*k1[2]*h))
zf = zi+vec(h*(c1*k1.+c2*k2));
zi = zf;
s = [s zf]
z1 = zi[1];
z2= zi[2];

end
s1 = s[1,:];
s2 = s[2,:];

plot(xi,s1)
plot!(xi,s2)



a = 0;
b = 2*pi;
n = 300;
h =(b-a)/n #length of the interval sections
xi = a:h:b #whole inerval of x

z1 = 1
z2 = 0
z3 = 0
z4 = 1
zi = [z1;z2;z3;z4]
s = zi;
f1(xi,z1,z2,z3,z4) = [z3;z4;-z1*(z1.^2+z2.^2).^(-3/2);-z2*(z1.^2+z2.^2).^(-3/2)]

for i in 1:n

k1 = f1.(xi[i],z1,z2,z3,z4)
k2 = f1.((xi[i]+h/2),(z1+k1[1]/2*h),(z2+k1[2]/2*h),(z3+k1[3]/2*h),(z4+k1[4]/2*h))
k3 = f1.((xi[i]+h/2),(z1+k2[1]/2*h),(z2+k2[2]/2*h),(z3+k2[3]/2*h),(z4+k2[4]/2*h))
k4 = f1.((xi[i]+h),(z1+k3[1]*h),(z2+k3[2]*h),(z3+k3[3]*h),(z4+k3[4]*h))
k = 1/6*(k1+2*k2+2*k3+k4)
zf = zi+k*h;
zi = zf;
s =[s zf]
z1 = zi[1];
z2 = zi[2];
z3 = zi[3];
z4 = zi[4];
end

s1 = s[1,:];
s2 = s[2,:];
s3 = s[3,:];
s4 = s[4,:];

using Plots
plotly()


plot(s1,s2)
scatter!(s3,s4)

scatter(s3,s4,aspect_ratio=:equal)



scatter(xi,s3)
scatter!(xi,s4)
