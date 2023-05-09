
#---------------linear splines----------------

x = [ 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0];
y = [ 0.0
  0.27182818284590454
  0.6847555777154057
  1.2769783442087028
  2.0935476878376944
  3.1874451224589215
  4.620817846279509
  6.466396377709597
  8.809119688943417
 11.747996543962481
 15.398235652779242];

# orginizing the data from xo -> xn in Ascending order
for i in 1:length(x)
    j = i;
    for i in 1:length(x)
if(x[j]>x[i])
    m = x[i];
    x[i] = x[j];
    x[j] = m;
    n = y[i];
    y[i]=y[j];
    y[j] = n;
end
end
end
x1 = x[length(x)]
y1 = y[length(y)]
m = length(x)-1;
for i in 1:length(x)-1
x1 = [ x1 x[m]];
y1 = [y1 y[m]];
m = m-1;
end
x = x1;
y = y1;
#calculating linear splines
f1 = 0;
using SymPy;
a = symbols("a")
for i in 1:length(x)-1
display("----------------")
println([ "for : " x[i] "<= x <=" x[i+1] ])
f3 = y[i]+(y[i+1]-y[i])*(x[i+1]-x[i]).^-1*(a-x[i])
display(f3)
f1 = [f1 f3]
end

f2 = f1[2]
for i in 3:length(f1)
f2 = [f2 f1[i]];

end

#------------------Plotting the linear spline--------------



f3 = [1.47220275232618*a - 0.0367031834250534 2.09701496538632*a - 0.224146847343094 2.92633129939635*a - 0.63880501434811]
f3 = vec(f3)
using Plots
plotly()
scatter!(x,y)


inter = 1:0.1:2
scatter(x,y)

for i in 1:length(x)-1

inter = x[i]:0.1:x[i+1]
f = f3[i]
plot!(inter,f.(inter),linewidth="3")

end
inter = x[length(x)-1]:0.1:x[length(x)]
plot!(inter,f.(inter),linewidth="3")

plot!(x,f.(),linewidth="3")

# f = a->f3[1]
plot!(inter,f.(inter),linewidth="3",color="red")

f = a->20

scatter!(x,y)
plot(x,y,linewidth="3",color="blue")



#-----------------Quadratic Splines ----------------------
using SymPy;

x = [0 10 15 20 22.5 30]
y = [0 227.04 362.78 517.35 602.97 901.67]

a = symbols("a"*string(1))
b = symbols("b"*string(1))
c = symbols("c"*string(1))

for i in 2:length(x)-1

a = [a symbols("a"*string(i))]
b = [b symbols("b"*string(i))]
c = [c symbols("c"*string(i))]

end
eq = a[1]*x[1].^2+b[1]*x[1]+c[1]-y[1];
for i in 2:length(x)-1

  eq = [eq a[i-1]*x[i].^2+b[i-1]*x[i]+c[i-1]-y[i]]
  eq = [eq a[i]*x[i].^2+b[i]*x[i]+c[i]-y[i]]
  if(i==length(x)-1)
      eq=[eq a[length(x)-1]*x[length(x)].^2+b[length(x)-1]*x[length(x)]+c[length(x)-1]-y[length(x)]]

  end
end
eq2 =0;
for i in 1:(length(x)-2)

eq2 = [eq2 2*a[i]*x[i+1]+b[i]-2*a[i+1]*x[i+1]-b[i+1]]

end

eq3 = eq2[2];

for i in 3:length(eq2)

eq3 = [eq3 eq2[i]]

end

eq3 = [eq3 a[1]]
eq = [eq eq3]

eq1 = eq[1]
eq2 = eq[2]
eq3 = eq[3]
eq4 = eq[4]
h =solve([eq[1],eq[2],eq[3],eq[4],eq[5],eq[6],eq[7],eq[8],eq[9],eq[10],eq[11],eq[12],eq[13],eq[14],eq[15]])

#--------------------------modified for matrix-----------------------------------------------

x = [0 10 15 20 22.5 30];
y = [0.0 227.04 362.78 517.35 602.97 901.67];
m = 1;
h = 2;
r = 1;
M = zeros((length(x)-1)*3,(length(x)-1)*3)

for i in 1:length(x)

b = [x[i]^2 x[i] 1]
if(i==1 || i==length(x))
    for g in 1:length(b)
        M[r,m] = b[g]
        m = m+1;
    end
    m = m-3;
    r = r+1;
end
if(i!=1 && i!=length(x))
b = [x[i]^2 x[i] 1]
for j in 1:length(b)
M[r,m] = b[j]
m = m+1;
end

r = r+1;
for j in 1:length(b)
M[r,m] = b[j]
m = m+1;
end
r = r+1;
m = m-3
end
end
m =1;
for i in 2:(length(x)-1)
b = [2*x[i] 1 0 -2*x[i] -1 0 ]

for j in 1:length(b)
M[r,m] = b[j]
m = m+1;
end

r = r+1;
m = m-3
end

M[(length(x)-1)*3,1] = 1;


y1 = y[1];
for h in 2:(length(y))

    for j in 1:2
        if(h!=length(y))
    y1 = [y1 y[h]]
    end
    if(h==length(y))
        y1 = [y1 y[h]]
        break;
    end
end
end

for g in 1:length(x)-1

y1 = [y1 0]

end
y1 = vec(y1);

a = M\y1



#-----------------------------Sparse matrix -------------------------------
using LinearAlgebra
using SparseArrays
using SymPy
x = [-2 -1.5 0 1.5]
y = [5 3 1 2]
a = symbols("a")

m = 1;
r =0;
for i in 1:(length(x)-1)*2
for j in 1:3
r = [r i]

end
end

c = 0;
m = 1;
h = 2;

for i in 1:(length(x)-1)*2
for j in 1:3
c = [c m]
m = m+1;
end
m = m-3;
if(i==h)
m = m+3;
h = h+2;
end
end

v = [x[1].^2 x[1] 1]
for i in 2:length(x)
if(i!=length(x))
for j in 1:2
v = [v x[i].^2 x[i] 1]
end
end
if(i==length(x))
v = [v x[i].^2 x[i] 1]
end
end
#converting columns
c1 = c[2]

for i in 3:length(c)

c1 = [c1 c[i]]
end

r1 = r[2]

for i in 3:length(c)

r1 = [r1 r[i]]
end

#s = sparse(r1,c1,v,15,15)


for i in (length(x)-1)*2+1:(length(x)-1)*3-1

for j in 1:6
r1 = [r1 i]
end


end
m = 1;
for i in 1:length(x)-2

for i in 1:6
c1 = [c1 m]
m = m+1
end
m = m-3;
end
c1 = [c1 1];
r1 = [r1 (length(x)-1)*3]

for i in 2:length(x)-1

v = [v 2*x[i] 1 0 -2*x[i] -1 0]

end
v = [v 1];
r1 = vec(r1)
c1 =vec(c1)
v = vec(v)
s = sparse(r1,c1,v,3*(length(x)-1),3*(length(x)-1))
A = Array(s)



y1 = y[1];
for h in 2:(length(y))

    for j in 1:2
        if(h!=length(y))
    y1 = [y1 y[h]]
    end
    if(h==length(y))
        y1 = [y1 y[h]]
        break;
    end
end
end

for g in 1:length(x)-1

y1 = [y1 0]

end
y1 = vec(y1);

Coef = A\y1

g = 0;

m = 1;
f1 = 0;
m1 = convert(Int64,length(Coef)/3)
for i in 1:length(x)-1

f = Coef[m]*a^2+Coef[m+1]*a+Coef[m+2]
f1 = [f1 f]
m = m+3;

end
f1 = vec(f1)
display(f1)




#---------------------------quadratic splines - recursion-----------------

x = [-1 0 0.5 1 2 2.5];
y = [2 1 0 1 2 3];
z = 0;
t = z;
for i in 1:length(x)-1

z =-z+2*(y[i+1]-y[i])/(x[i+1]-x[i])
t = [t z]
end

using SymPy
z = t;
a = symbols("a")
p = 0;
for i in 1:length(z)-1

q = (z[i+1]-z[i])*1/2*(x[i+1]-x[i]).^-1*(a-x[i]).^2+z[i]*(a-x[i])+y[i]
p = [p q]
end
p = vec(p);
display(p)

#-------------------------Cubic splines-----------------------------
#---------rows------------
using LinearAlgebra
using SparseArrays
using SymPy

x = [1 2 3 4 5]
y = [0 1 0 1 0]


r = 0;

if(length(x)>=3)
val = 3;
k = 1;
for i in 1:100
if(length(x)!=val)
val = val+1;
k = k+3;
end
if(length(x)==val)
break;
end
end
end

for i in 1:k
if(i!=1 && i!=k)

for j in 1:3

    r = [r i]

end
end
if(i==1 || i==k)
for m in 1:2
r = [r i]
end
end
end
r1 = r[2]
for i in 3:length(r)
r1 = [r1 r[i]]

end
r2 = 0;

for i in 1:k

r2 = [r2 r1[i]]

end
r3 = r2[2]
for i in 3:length(r2)
r3 = [r3 r2[i]]

end

#r3 final



# ----------columns-------------
c = 0;
m = 1;
n = 0;
for i in 1:k
if(i==1 || i==k)
    for j in m:m+1
        c = [c j]
end
end
m = length(x)-3
if(i!=1 && i!=k)
    for g in n:n+2

        c = [c g]

    end

end
n = n+1;
end

c1 = c[2];

for i in 3:length(c)
c1 = [c1 c[i]]

end
c2 = 0;
for i in 1:k
c2 = [c2 c1[i]]

end

c3 = c2[2]

for i in 3:length(c2)

c3 = [c3 c2[i]]
end

#c3 final


#---------values -----------
h1 = 0;
for i in 1:length(x)-1
h1 =[h1 x[i+1]-x[i]]

end
h = h1[2]
for i in 3:length(h1)
h = [h h1[i]]

end
h1 = h;

#h1 final
#-------------------
u  =0;
for i in 1:length(h)-1
u =[u 2*(h[i]+h[i+1])]

end
u1 = u[2]
for i in 3:length(u)

u1 = [u1 u[i]]

end

#u1 final

#-----------------------------------------------------------------
v = 0;
for i in 1:length(u1)

    v = [v u1[i]]
    if(i==length(u1))
        break;
    end
    for j in 1:2
        v = [v h[i+1]]

    end

end


v1 = v[2]

for i in 3:length(v)
v1 = [v1 v[i]]
end

#final v1


if(length(v1)>1)
v1 = vec(v1);
r3 = vec(r3);
c3 = vec(c3);
end
if(length(v1)==1)
    v1 = [v1];
    r3 = [r3];
    c3 = [c3];
end
s = sparse(r3,c3,v1,(length(x)-2),(length(x)-2))
A = Array(s)

b = 0;
for i in 1:length(h)

b = [b (y[i+1]-y[i])/(h[i])]

end
b1 = b[2]

for i in 3:length(b)
b1 = [b1 b[i]]
end

m = 0;

for i in 1:length(b1)-1

m =[m 6*[b1[i+1]-b1[i]]]

end
m1 = m[2]
if(length(m)>=3)
for i in 3:length(m)
m1 = [m1 m[i]]
end
end
if(length(m1)>1)
m1 = vec(m1)
end
if(length(m1)==1)
m1 = [m1];
end
#solution for z ->
z = A\m1

a = symbols("a")
t = 0;
z1 = 0
for i in 1:(length(z))
z1 =[z1 z[i]]
if(i==length(z))
z1 = [z1 0]
end
end
z = z1;
for i in 1:length(z)-1

p = z[i+1]/(6*h[i])*(a-x[i]).^3+z[i]/(6*h[i])*(x[i+1]-a).^3+(y[i+1]/h[i]-h[i]*z[i+1]/6)*(a-x[i])+((y[i]/h[i])-(h[i]/6*z[i]))*(x[i+1]-a)
t = [t p];
end

t1 = t[2]

for i in 3:length(t)
t1 =  [t1 t[i]]
end

t1 = vec(t1)



#plotting

scatter(x,y)
inter = 0.3:0.01:0.4
plot!(inter, f3.(inter),linewidth="3")
x1 = [0.18 0.22 0.34]
y1 = [-0.3487118399999998 -0.221380 0.107142]

scatter!(x1,y1)

y = [0.86199480 0.95802009 1.0986123 1.2943767]
x = [-1 -0.5 0 0.5]
