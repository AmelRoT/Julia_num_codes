##problem 1 a
x = 0.2:0.0001:0.3
y = ((x.*(cos.(x)))-2*(x.^2)+3*x).-1.0000

using Plots
plotly()
plot(x,y,linewidth=3,color="blue")

## problem 2c

x = 0:0.0001:5;
y = (log.(x)-x.^2+(5/2)*x).-1.0;
y = tand.(x)
plot(x,y,linewidth=3)

## problem 5

x = -20:0.0001:20
y = 2*x.*cos.(2*x)-(x.-2).^2
plot(x,y,linewidth=3)

x1 = 0:0.1:1;
y1 = (-exp.(x1)+((exp(1)-1)*sin.((pi/2)*x1))).+1
counter = 0;
for i in 1 : length(x1)

if (round(y1[i],digits=3) == 0.0)
    counter  = counter+1;
    display(counter)
end

end
#-------------------------------------------------------------------
x = 0:0.001:1
y = ((x.-1.0).*tan.(x)).+x.*sin.(pi*x)

plot(x,y,linewidth=3)

#--------------------------------------------------------------
x = -5:0.0001:5
y = (x.^3+4x.^2).-10
y1=(x1.^3+4x1.^2).-10
plot(x,y,linewidth=3)

bsc = 1;
uplimit = 2
lowerlimit = 0
x = uplimit;
x1 = lowerlimit;

for i in 1:10
    bsc = (uplimit-lowerlimit)/2

if(y>0 && y1 <0)
    lowerlimit = bsc;
end

end



#------------------------------------------------------------------------------
x1 = 2;
x = 0;
y = (x.^3+4*x.^2).-10;
y1=(x1.^3+4*x1.^2).-10;
b = 0;
for i in 1:300
    b = (x1-x)/2
    x = x+b;
    y = (x.^3+4x.^2).-10;
if (y < 0 )
    display(x)
end
if(y>0)
x1 = x
x = 1;

end

if(y==0)
    break;
end
end

# -------------------------------------------------------------------------
x = -2:0.0001:2
f = x->(x.^5+3*cos.(x)+exp(2*log.(x))+2*atan.(x)).-10
using Plots
plotly()
plot(x,y, linewidth=3,color="blue")
#------------------------not in usage now --------------------------------------
f_a = (8*a.^3+4*a.^2+4*a).-15
f_b = (8*b.^3+4*b.^2+4*b).-15

for i in 1:length(a)
 if((f_a[i])*(f_b[i])<0)
     a = a[i]
     b = b[i]
     f_a = f_a[i]
     f_b = f_b[i]
     break;
 end
end
#----------------------------------------------------------------------------------
# modified version of the given code
# ----------------------------Bisection Method -----------------------------------
f = x->(-6.8*x.^2+32.17*(2*sin.(x)-exp.(x)+exp.(-x)))
#--------a--------------
f = x-> x^3-25
#-------b---------------
f = x-> x^2-3
f = h->12.4-10*((0.5*pi-asin(h))-h*(1-h.^2)^(1/2))
f = x-> x.^3+4x.^2-10
a = -1:0.0001:0
b =  0:0.0001:1


function BisectionMethod(f,a,b)

for i in 1:length(a)
 if(sign(f.(a[i]))*sign(f.(b[i]))<0)
     a = a[i]
     b = b[i]
     break;
 end
end
    f_a = f.(a)
    f_b = f.(b)

# f(a)*f(b)<0 condition


# the number of iterations ->

epsi = 10^-10

    n = (log(abs(a-b)/epsi)).*log(2)^-1
    n = round(n,digits=0)

#               finding the roots
for i in 1:n

bis = a+((b-a)/2);
f_bis = f.(bis);
#                    First case scenario
    if(sign(f_bis)*sign(f_a) > 0 && sign(f_a) < 0 )
        a = bis;
        f_a = f.(a);
        display(a)

    end

#                    Second case scenario
        if (sign(f_bis)*sign(f_a) > 0 && sign(f_a) > 0 )
        a = bis;
        f_a = f.(a);
        display(a)
    end

#                    Third case scenario
    if(sign(f_bis)*sign(f_b) > 0 && sign(f_b) < 0 )
        b = bis;
        f_b =f.(b);
        display(b)
    end
    #                Fourth case scenario
    if(sign(f_bis)*sign(f_b) > 0 && sign(f_b) > 0 )
        b = bis;
        f_b = f.(b);
        display(b)
    end
                #    Fifth case scenario
    if(f_bis==0)
    root = bis;

    end

end
println("Number of iterations : ", n)
end

#----------------------------------- End ------------------------------------------
#-----------c-------------------------------
using Plots
f =x-> exp(x)-2 # result -> 0.69314718....
f = x->cos.(exp.(x)-2)
# for the cos interval ->1.27272788...
f = x-> cos.(exp(x)-2)+2-exp.(x) #1.00762397155....
a = 0.5
b = 1.5
x1 = 0.5:0.0001:1.5
plotly()
plot(x1,f,linewidth=3)

#--------problem 6--------------------------------------
using SymEngine
L = symbols(:L)
r = symbols(:r)
h = symbols(:h)
theta = symbols(:theta)
V = L*(0.5*pi*r.^2-r.^2*theta-h*(r.^2-h.^2).^(1/2))
L = 10
r = 1
V = 12.4
h = 0.01
f = h->12.4-10*((0.5*pi-asinx(h))-h*(1-h.^2)^(1/2))
f = h->sin(0.33079642-h*(1-h.^2)^(1/2))-h
x=-1:0.0001:1
plot(x,f,linewidth=3)
root_h=BisectionMethod.(f,-1,1)
depth = 1-root_h
#find the h -> depth = r-h  = 1-0.1652234375... = 0.8347765625...

#----------------------Problem 7-------------------------
f = w-> -1.7-(1/2*(32.17)*w.^-2*(((exp.(w*1)-exp.(-w*1))/2-sin.(w*1))))
t = 1
plot(x,f,linewidth=3)
BisectionMethod.(f,-0.01,-1)

f =x-> x^5+2*x^3-5x-2
x = 0:0.0001:2.5
using Plots
plot(x,f,linewidth=3,color="blue")
plotly()



#----------------------------------PP----------------------------
#9->
y = x->exp(x)-2
y1 = x-> cos(exp(x)-2)
x = 0:0.0001:2.5
plot!(x,y,linewidth=3)
plot!(x,y1,linewidth=3)

f = x-> exp(x)-2
BisectionMethod(f,0,3)

g = x-> exp(x)-2-cos(exp(x)-2)
BisectionMethod(g,0,4)

#next one
f = x->(x+2)*(x+1)*x*(x-1)^3*(x-2)

a = -2.5
b = 3
BisectionMethod(f,a,b)
#applied one
f = h-> sin(-1.24+1/2*pi-h*(1-h^2).^(1/2)-h)
a = 0
b = 1
BisectionMethod(f,a,b)

f1 = h->-12.4+10*(1/2*pi-asin(h)-h*(1-h^2)^(1/2))


#another applied one
f = w->-1.7-(32.17/2)*(w^-2)*(exp(w)-exp(-w)-2*sin(w))/2
plot(x,f,linewidth=3)
BisectionMethod(f,-2,1)
x = -2:0.0001:2

#---------------------------PP Newt-----------------------------
# newton method fast version
using Calculus
using ForwardDiff
using ChainRules
xo = 1;
f= x-> x-10
f1 = derivative(f)
function Newton(f,f1,xo)

for i in 1:200
x = xo-(f.(xo)*((f1.(xo))).^-1)
if (abs(x-xo)<=10^-10)

display(x)
display(i)
end
xo = x

end

end

#rapid secant method

function Secanta(f,x0,x2)

for i in 1:200

x1 = x0-((f.(x0)*(x0-x2))*(f.(x0)-f.(x2))^-1)
x2 = x0;
if(abs(x1-x0) <10^-10)
display(x1)
display(i)
break;

end
x0 = x1;


end





end

f = x ->2-exp.(x);
xo = 3;
x1= 2;


##
using Plots;
plotly();
f = x-> ((x-1)*(x+1))/(x+2)
x1 = -10:0.05:10;
plot(x1,f.(x1),linewidth=3)

#

f = x-> ((x-2)*(x-3)*x^2)/(x+1)
x1 = -3.5:0.05:10;
plot(x1,f.(x1),linewidth=3)


f = x-> ((x+4)*(2x-1))/(x.^2-2*x-3)
x1 = -6:0.05:10;
plot(x1,f.(x1),linewidth=3)


#-------------------------------------------------#
A = 858/(sind(38))
B = 1174.89/(sind(58))
C=  1378.25/sind(84)

F = (200*cosd(30))/cosd(45)
Fr = F*sind(45)+200*sind(30)

Fbx = 2*cosd(33.69)
Fby = 2*sind(33.69)
Fb = Fbx^2+Fby^2

Frx = Fbx+3;
Fry = Fby
FR = sqrt(Frx^2+Fry^2)

g = cosd(30)/sind(30)

Fr = sind(78.6)/(sind(30))*2

atand(3/4)

9.8*6

T21 = (cosd(45)/cosd(36.869)*sind(36.869)+sind(45)).^-1*58.8
T3 = T21 *cosd(45)/cosd(36.869)

T3 = 80/sind(30)

atand(4/3)

(100*sind(53.13)+100*sind(30))/9.8
a = 1/sqrt(16+2*144)

b = [4 12 -12]
c = a*b*2

d1 = (sqrt(2)+sqrt(6))/(4)

t1 = atand(0.2/0.4)
t2 = 90-t1

Ay = (+225*0.1-120*0.125)/0.4

Mc = 18.75*0.25+40*0.025

sqrt(2^2+1.5^2)


Pa = 87*10^3+8720*6*10^-2+7*133100*10^-2-9790*10^-2*5
Ixx = 1/4*pi*25^4+25^2*pi*75^2

IiX = 150^3*100/12+(100*(75*2)*75^2)

Moment = IiX-Ixx

Fh = 62.4*12*24*50

z1 = ((24)^(5/2))/((24)^(3/2))*3/5

x1 = 575/(156.66666)
