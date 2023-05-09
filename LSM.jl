#---------------------------First order Diff with 2 points----------------------------
using Plots

x = 4.0:0.2:8.0
y = [-5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59];
d1 = 0;
function dx(x,y)
d1 = 0;
d = 0;

for i in 1:length(x)


if(i==1)
d = [d (y[i+1]-y[i])/(x[i+1]-x[i])]

end


if(i!=1 && i!=length(x))

d = [d (y[i+1]-y[i-1])/(x[i+1]-x[i-1])]

end



if(i==length(x))
d = [d (y[i]-y[i-1])/(x[i]-x[i-1])]
end

end

d1 = d[2]

for i in 3:length(d)

d1 = [d1 d[i]]
end

d1 = vec(d1)
return d1
end
d1 = dx(x,y)
d2 = dx(x,d1)
scatter(x,y)
scatter!(x,d1)
scatter!(x,d2)

#---------------------------First order Diff with 3 points----------------------------
using Plots

x = [1 2 3 4 5]
y = [3 4 5 6 7];

x = [7.4 7.6 7.8 8]
y = [-68.3193 -71.6982 -75.1576 -78.6974]
f =x-> log.(x+2)-(x+1).^2
f2 = derivative(f)

d = 0;
for i in 1:length(x)


if(i==1)

d = [d (-3*y[i]+4*y[i+1]-y[i+2])/(2*(x[i+1]-x[i]))]


end



if(i==2 || i==length(x)-1)

d = [d (y[i+1]-y[i-1])/(x[i+1]-x[i-1])]
#d = [d (y[i+1]-y[i-1])/(x[i+1]-x[i-1])]

end

if(i==length(x))

d = [d (y[i-2]-4*y[i-1]+3*y[i])/(2*(x[i]-x[i-1]))]

end



if(i!=1 && i!=length(x) && i!=2 && i!=length(x)-1)
d = [d (y[i+1]-y[i-1])/(x[i+1]-x[i-1])]
#d = [d (y[i-2]-8*y[i-1]+8*y[i+1]-y[i+2])/(12*(x[i+1]-x[i]))]
end

end

d1 = d[2]

for i in 3:length(d)

d1 = [d1 d[i]]
end

d1 = vec(d1)
scatter!(x,y)
scatter!(x,d1)

j = x->2*exp.(2*x)+2*sin(2*x)
k = derivative(j)

for i in 1:length(x)

k1 = abs(d1[i]-j(x[i]))
display(k1)
end


#---------------------------Second order Diff with 2 points----------------------------
using Plots

x = [1 2 3 4 5]
y = [3 4 5 6 7];
d = 0;

for i in 1:length(x)


if(i==1)

d = [d (y[i]-2*y[i+1]+y[i+2])/(x[i+1]-x[i]).^2]


end




if(i==length(x))

d = [d (y[i-2]-2*y[i-1]+y[i])/(x[i]-x[i-1]).^2]

end



if(i!=1 && i!=length(x))
d = [d (y[i-1]-2*y[i]+y[i+1])/(x[i+1]-x[i]).^2]
end

end

d1 = d[2]

for i in 3:length(d)

d1 = [d1 d[i]]
end

d1 = vec(d1)
scatter!(x,y)
scatter!(x,d1)


#---------------------------Second order Diff with 3 points----------------------------
using Plots

x = [1 2 3 4 5]
y = [3 4 5 6 7];

x = [1.2 1.29 1.3 1.31 1.4]
y = [11.59006 13.78176 14.04276 14.30741 16.86187]

f =x-> (3*x*exp.(x)).-cos.(x)
f2 = derivative(f)
f3 = derivative(f2)
d = 0;

for i in 1:length(x)


if(i==1)

d = [d (2*y[i]-5*y[i+1]+4*y[i+2]-y[i+3])/(x[i+1]-x[i]).^2]


end



if(i==2 || i==length(x)-1 ||i==length(x)-2)

d = [d (y[i-1]-2*y[i]+y[i+1])/(x[i+1]-x[i]).^2]

end


if(i==length(x))

d = [d (-y[i-3]+4*y[i-2]-5*y[i-1]+2*y[i])/(x[i]-x[i-1]).^2]

end



if(2<i<(length(x)-2))
d = [d (-y[i-2]+16*y[i-1]-30*y[i]+16*y[i+1]-y[i+2])/(12*(x[i+1]-x[i]).^2)]
end

end

d1 = d[2]

for i in 3:length(d)

d1 = [d1 d[i]]
end

d1 = vec(d1)
scatter!(x,y)
scatter!(x,d1)



#-----------------------Just Difference---------------

h = [0.1 0.01 0.001]
xi = 1
g = x-> log.(x)

#quiz 3
#0.9531017980432493
#0.9950330853168092
#0.9995003330834232

#forward difference
for i in 1:length(h)
f1 = (g.(xi+h[i])-g.(xi))/h[i]
display(f1)
end


#backward

for i in 1:length(h)
f1 = (g(xi)-g.(xi-h[i]))/h[i]
display(f1)
end


#central

for i in 1:length(h)
f1 = (g.(xi+h[i])-g.(xi-h[i]))/(2*h[i])
display(f1)
end

#--------------Using formulas for first order---------------------------------

x = [6.0 6.1 6.2 6.3 6.4]
y =[0.1750 -0.1998 -0.2223 -0.2422 -0.2596]

x = [0.5 0.6 0.7]
y = [0.4794 0.5646 0.6442]

x = [0 0.2 0.4]
y = [0 0.74140 1.3718]

x = [-0.3 -0.2 -0.1]
y = [1.9507 2.0421 2.0601]

x = [1 1.2 1.4]
y = [1 1.2625 1.6595]



f = x->x.^2*log.(x)+1
f2 = derivative(f)
#forward
display(" ---- ")
for i in 1:length(x)-1
f = (y[i+1]-y[i])/(x[i+1]-x[i])
display(f)

end

display(" ---- ")
#backward
for i in 2:length(x)
f = (y[i]-y[i-1])/(x[i]-x[i-1])
display(f)
end


display(" ---- ")
#central


for i in 2:length(x)-1

f = (y[i+1]-y[i-1])/(x[i+1]-x[i-1])
display(f)
end

#-------------------------Using formulas for second order--------------------------------------

x = [0.5 0.6 0.7]
y = [0.4794 0.5646 0.6442]

x = [1 2 3 4 5]
y = [2.4142 2.6734 2.8974 3.0976 3.2804]
d = 0

for i in 3:3
d = [d (y[i-2]-8*y[i-1]+8*y[i+1]-y[i+2])/(12*(x[i+1]-x[i]))]
end


x = [0.2 0.4 0.6 0.8 1.0]
y = [0.9798652 0.9177710 0.808038 0.6386093 0.3843735]

h = 0.001
h = 0.1
xi = 1.0
f =x->(x).^-1
d = (f.(xi-h)-2*f.(xi)+f.(xi+h))/(h.^2)

#2.0202020202020217
#2.0002000200025627
#2.0000020002353125


#-----------------------------Least Square Approximation Linear------------------------------------------

x = [0 10 20 30 40 50 60 70 80 90 100]
y = [0.94 0.96 1.0 1.05 1.07 1.09 1.14 1.17 1.21 1.24 1.28]

sumX = sum(x);
sumY = sum(y);

xy = x[1]*y[1];
for i in 2:length(x)

xy = [xy x[i]*y[i]];

end
xy1 = sum(xy)
xsq_1 = x.^2
xsq = sum(xsq_1)
m = length(x)

M = [sumX m;xsq sumX]

R = [sumY;xy1]

a =M\R

using SymPy

s = symbols("s");

y1 = a[1]*s+a[2]

#--------------------------------Quadratic Square Approximation ---------------------------
using SymPy
x = [1.2 1.4 1.6 1.8 2.0]
y = [0.504545 0.524791 0.557778 0.602524 0.659112]

sumX = sum(x);
sumY = sum(y);

xy = x[1]*y[1]
for i in 2:length(x)

xy = [xy x[i]*y[i]];

end
xy1 = sum(xy)
xsq_1 = x.^2
xsq = sum(xsq_1)
xcub_1 = x.^3
xcub = sum(xcub_1)
xfour_1= x.^4
xfour = sum(xfour_1)
xsqy_1 = xsq_1[1]*y[1]

for i in 2:length(x)

xsqy_1 = [xsqy_1 xsq_1[i]*y[i]];

end
xsqy = sum(xsqy_1)

m = length(x)

M = [m sumX xsq;sumX xsq xcub;xsq xcub xfour]

R = [sumY;xy1;xsqy]

a = M\R

s = symbols("s")

y1 = a[1]+a[2]*s+a[3]*s^2

#----------------------------------------------------------------------------
for i in 1:length(x)

g = [g (abs(y1(i)-y[i]))]


end
a = [
  0.230769230769226
  11.4615384615385
  17.3076923076923
  15.9230769230769
  2.15384615384616
  9.38461538461539
  12.3846153846154
  25.1538461538462
  4.07692307692308
  7.69230769230769
  5.53846153846154
  13.7692307692308]

x = [1 2 3 4 5 6 7 8 9 10 11 12]
y = [78 65 92 57 69 60 80 91 60 70 55 45]
f = a->80.0 - 1.76923076923077*a
for i in 1:length(x)

g = [g (abs(f(i)-y[i])).^2]


end

x = [0 2 4 6 9 11 12 15 17 19]
y = [5 6 7 6 9 8 7 10 12 12]
f = s->0.35246*s + 4.85153

for i in 1:length(x)

g = [g ((-f(x[i])+y[i])).^2]


end

display(sum(g))
