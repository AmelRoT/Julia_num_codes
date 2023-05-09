#-------------------Trapezoid Rule--------------------


x = [1.8 2.0 2.2 2.4 2.6]
y = [3.12014 4.42569 6.04241 8.03014 10.46675]
Iby = (2.6-1.8)/4;
x = 1.8:Iby:2.6
#interval ->[a,b]

a = 0;
b = 2;

h = b-a
#number of points (dividing the interval)
n = 100;

#stepping the interval by

a = 0
b = 3
Iby = (b-a)/n

x = a:Iby:b

y = x.^2

T = 0;
for i in 1:length(x)

if(x[i]==a)

T =[T Iby/2*(y[i])]

end
if(x[i]==b)

T =[T Iby/2*(y[i])]

end

if(x[i]!=a && x[i]!=b)

T =[T Iby/2*(2*y[i])]

end

end

T1 = T[2]

for i in 3:length(T)

T1 = [T1 T[i]]

end

A = sum(T1)



#------------------------------------------
#Trapezoid method --- remade version for easier writing


x = [1.8 2.0 2.2 2.4 2.6]
y = [3.12014 4.42569 6.04241 8.03014 10.46675]
Iby = (2.6-1.8)/4;
x = 1.8:Iby:2.6
#interval ->[a,b]

a = 0;
b = 2*pi;

h = b-a
#number of points (dividing the interval)
n = 500;

#stepping the interval by
Iby = h/n

x = a:Iby:b
y = x-> (sin(7*x)).^2

T = 0;
for i in 1:length(x)

if(x[i]==a)

T =[T Iby/2*(y.(x[i]))]

end
if(x[i]==b)

T =[T Iby/2*(y.(x[i]))]

end

if(x[i]!=a && x[i]!=b)

T =[T Iby/2*(2*(y.(x[i])))]

end

end

T1 = T[2]

for i in 3:length(T)

T1 = [T1 T[i]]

end

A = sum(T1)




#----------------------Simpson's rule 1/3--------------------------


x = [0 6 12 18 24 30 36 42 48 54 60 66 72 78 84]
y = [124 134 148 156 147 133 121 109 99 85 78 89 104 116 123]

#interval ->[a,b]
Iby = 0.05

a = 0;
b = 24;

h = b-a
# number of points (dividing the interval)

n = 500;

#stepping the interval by
Iby = h/n

x = a:Iby:b
y = x->  2.041x^(1/2)

T = 0;

for i in 1:length(x)

if(x[i]==a)

T = [T Iby/3*y.(x[i])]


end

if(x[i]==b)

T = [T Iby/3*(y(x[i]))]

end
if(x[i]!=b && i%2==0)

T = [T Iby/3*(4*y(x[i]))]

end

if(x[i]!=a && x[i]!=b && i%2==1)

T = [T Iby/3*(2*y(x[i]))]

end

end

T1 = T[2]

for i in 3:length(T)

T1 = [T1 T[i]]


end

A = sum(T1)



#----------------------Simpson's rule 1/3--------E-Writing----------(use it to solve definite integrals)

x = [0 6 12 18 24 30 36 42 48 54 60 66 72 78 84]
y = [124 134 148 156 147 133 121 109 99 85 78 89 104 116 123]


#interval ->[a,b]
Iby = 0.05

a = 0;
b = 2;

h = b-a
# number of points (dividing the interval)

n = 100;

#stepping the interval by
Iby = h/n

x = a:Iby:b
y =x->x^2

T = 0;

for i in 1:length(x)

if(x[i]==a)

T = [T Iby/3*(y.(x[i]))]

end


if(x[i]==b)

T = [T Iby/3*(y.(x[i]))]

end
if(x[i]!=b && i%2==0)

T = [T Iby/3*(4*y.(x[i]))]

end

if(x[i]!=a && x[i]!=b && i%2==1)

T = [T Iby/3*(2*y.(x[i]))]

end

end

T1 = T[2]

for i in 3:length(T)

T1 = [T1 T[i]]

end

A = sum(T1)





#------------------------------Problems------------------------------

y = 200*(x.*(x.+7).^-1)*exp.(-2.5*x/30)



(0.27*10^3*9.81)/(13550*9.81)









#----------------------Simpson's rule 3/8--------------------------

x = [0 0.5 1 1.5 2 2.5 3]
y = [1 0.667 0.5 0.4 0.3333 0.2857 0.25 ]
Iby = 0.5
#interval ->[a,b]
a = 0;
b = pi/2;

h = b-a
# number of points (dividing the interval)

n = 3;

#stepping the interval by
Iby = h/n

x = a:Iby:b
y =exp.(sin.(x))

T = 0;
m = 2;
m2=4
for i in 1:length(x)

if(x[i]==a)

T = [T 3*Iby/8*(y[i])]

end

if(x[i]==b)

T = [T 3*Iby/8*(y[i])]

end
if(x[i]!=a && x[i]!=b && (i==m))

T = [T 3*Iby/8*(3*(y[i]+y[i+1]))]
m = m+3;
end

if(x[i]!=a && x[i]!=b && i==m2)

T = [T 3*Iby/8*(2*y[i])]
m2 = m2+3;
end

end

T1 = T[2]

for i in 3:length(T)

T1 = [T1 T[i]]


end

A = sum(T1)

















#----------------------Rectangle rule (Midpoint)--------------------------



#interval ->[a,b]

a = 0;
b = 3;

h = b-a
# number of points (dividing the interval)

n = 15 ;

#stepping the interval by
Iby = h/n

x = a:Iby:b
x1 = 0
for i in 1:n

x1 = [x1 (x[i]+x[i+1])/2]

end

x2 = x1[2];
for i in 3:length(x1)
x2 = [x2 x1[i]]

end

#use x2 for values in y

y =x->
T = 0;

for i in 1:length(x2)

T = [T Iby*y.(x2[i])]

end


T1 = T[2]

for i in 3:length(T)

T1 = [T1 T[i]]


end

A = sum(T1)

#---------------------------New one ---------------------------

a = 0
b = pi/2
n = 16
h = (b-a)/n
xx = zeros(1,n);
for i in 1:n
xx[i] = a+h/2+(i-1)*h
end
sum1= 0
for i in 1:n
sum1 = sum1+h*(exp(xx[i])-1)/sin.(xx[i])
end
display(sum1)







#----------------Gaussian Quadrature --------------

a = 0;
b = 3;

f =x->sqrt(((5*x.^2-81)/(x.^2-9)))/3

#------------changing the interval------------

dx = 1/2*(b-a)
t =x-> ((b-a)*x+a+b)*1/2

f1 =x-> dx*f.((((b-a)*x+a+b)*1/2))
f2  = x->dx*f.(t.(x))


#for 1.st n = 2

f3 = f1.(-1/sqrt(3))+f1.(sqrt(3)/3)


#for the 2nd n = 3
C1 = 0.5555556
C2 = 0.8888889
C3 = C1;

f3 = C1*f1.(-0.77459667)+C2*f1.(0)+C3*f1(+0.77459667)



#For the 3rd n = 4
C1 = 0.3478548
C2 = 0.6521452
C3 =C2;
C4 = C1;

f3 = C1*f1.(-0.86113631)+C2*f1.(-0.33998104)+C3*f1.(0.33998104)+C4*f1.(0.86113631)
#using t inside the function
f4 = C1*f2.(-0.86113631)+C2*f2.(-0.33998104)+C3*f2.(0.33998104)+C4*f2.(0.86113631)




#for the n = 5
C1 = 0.2369268851
C2 = 0.4786286705
C3 =C2;
C4 = C1;
C5 = 0.5688888889

f3 = C1*f1.(-0.9061798459)+C2*f1.(-0.5384693101)+C3*f1.(0.5384693101)+C4*f1.(0.9061798459)+C5*f1.(0)






#------------------------Trapezoid method for double integrals------------------------------
using Plots
plotly()
using SymPy

f(x,y) = x+y.^2

a = 1
b  = 2
c = 2
d = 3

n = 100
m = 100

k = (d-c)/n
h = (b-a)/m

g =(h/2)*(k/2)

t1 = a:h:b
t2 =c:k:d

V = 0;
for i in 1:length(t1)

for j in 1:length(t2)

if((t1[i]==a || t1[i]==b) && (t2[j]==c || t2[j]==d))
V = [V f(t1[i],t2[j])]

end

if((t1[i]!=a && t1[i]!=b) && (t2[j]!=c && t2[j]!=d))
V = [V 4*f(t1[i],t2[j])]

end

if((t1[i]==a || t1[i]==b) && (t2[j]!=c && t2[j]!=d))
V = [V 2*f(t1[i],t2[j])]

end

if((t1[i]!=a && t1[i]!=b) && (t2[j]==c || t2[j]==d))
V = [V 2*f(t1[i],t2[j])]

end

end
end

V1 = sum(V)*g



#------------------------1/3 Simpsons method for double integrals------------------------------
using Plots
plotly()
using SymPy

f(x,y) =x+y.^2
#x bounds
a = 1
b  = 2
#y bounds
c = 2
d = 3

n = 50
m = 50

k = (d-c)/n
h = (b-a)/m

g =(h/3)*(k/3)

t1 = a:h:b
t2 =c:k:d

V = 0;
for i in 1:length(t1)

for j in 1:length(t2)

if((t1[i]==a || t1[i]==b) && (t2[j]==c || t2[j]==d))
V = [V f(t1[i],t2[j])]

end

if((t1[i]==a || t1[i]==b) && (t2[j]!=c && t2[j]!=d))

if(j%2==0)
V = [V 4*f(t1[i],t2[j])]
end
if(j%2==1)
V = [V 2*f(t1[i],t2[j])]
end
end

#
if((t1[i]!=a && t1[i]!=b) && (t2[j]!=c && t2[j]!=d) && i%2==0)
if(j%2==0)
V = [V 16*f(t1[i],t2[j])]
end
if(j%2==1)
V = [V 8*f(t1[i],t2[j])]
end
end






if((t1[i]!=a && t1[i]!=b) && (t2[j]!=c && t2[j]!=d) && i%2==1)
if(j%2==0)
V = [V 8*f(t1[i],t2[j])]
end
if(j%2==1)
V = [V 4*f(t1[i],t2[j])]
end
end

if((t1[i]!=a && t1[i]!=b) && (t2[j]==c || t2[j]==d) && i%2==0)
V = [V 4*f(t1[i],t2[j])]

end
if((t1[i]!=a && t1[i]!=b) && (t2[j]==c || t2[j]==d) && i%2==1)
V = [V 2*f(t1[i],t2[j])]

end

end
end

V1 = sum(V)*g







#--------remake 1----------------1/3 Simpsons method for double integrals------------------------------
using Plots
plotly()
using SymPy

f(x,y) =x+y.^2
#x bounds
a = 1
b  = 2
#y bounds
c =x->x
d =x->x.^2

n = 4
m = n

k =x->((x.^2-x)/n+epsi)/n #(d-c)/n
h = (b-a)/m

g =(h/3)
epsi = 0.00000001;
t1 = a:h:b
t2 = x-> x:((x.^2-x)/n+epsi):x.^2

V = 0;
for i in 1:length(t1)

if((t1[i]==a || t1[i]==b))
t2_p = t2.(t1[i])

for j in 1:(n+1)
if(j==1 || j==(n+1))
V = [V (k.(t1[i])/3)*f.(t1[i],t2_p[j])]

end
if(j%2==0 && j!=1 && j!=(n+1))
V = [V (4*k.(t1[i])/3)*f.(t1[i],t2_p[j])]

end

if(j%2==1 && j!=1 && j!=(n+1))
V = [V (2*k.(t1[i])/3)*f.(t1[i],t2_p[j])]
end
end
end


if((t1[i]!=a && t1[i]!=b && i%2==0))
t2_p = t2.(round(t1[i],digits=5))
display(i)
for j in 1:(n+1)
if(j==1 || j==(n+1))

V = [V (4*k.(t1[i])/3)*f.(t1[i],t2_p[j])]
end
if(j%2==0 && j!=1 && j!=(n+1))
V = [V (16*k.(t1[i])/3)*f.(t1[i],t2_p[j])]
end

if(j%2==1 && j!=1 && j!=(n+1))
V = [V ((8*k.(t1[i])/3)*f.(t1[i],t2_p[j]))]
end
end
end


if((t1[i]!=a && t1[i]!=b && i%2==1))
t2_p = t2.(round(t1[i],digits=5))
for j in 1:(n+1)
if(j==1 || j==(n+1))
V = [V ((2*k.(t1[i])/3)*f.(t1[i],t2_p[j]))]
end
if(j%2==0 && j!=1 && j!=(n+1))
V = [V (8*k.(t1[i])/3)*f.(t1[i],t2_p[j])]
end

if(j%2==1 && j!=1 && j!=(n+1))
V = [V (4*k.(t1[i])/3)*f.(t1[i],t2_p[j])]
end
end
end
end
V1 = sum(V)*g
