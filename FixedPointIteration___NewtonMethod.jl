# Fixed point iteration method
using Plots

# --------------------Plotting a function -----------------------------------------------

plotly()
x = -2:0.0001:2
y = (x.^3+x).-1
y1 = x
plot(x,y,linewidth=3,color="blue")
plot!(x,y1,linewidth=3,color="red")
#-------------------------------------------------
f = x -> (x.^3+4*x.^2).-10 # any function that a user wants
f1 = x-> x # standard y = x function
test = 1.449
xo  = 1.25; # starting point or the fixed point
g =x-> ((x-2*x.^2).+3)^(1/4)

# --------------------------this is just for testing ---------------------------
for i in 1:length(test)

    if(f.(test[i]) == f1.(test[i]))
        xo = test[i]
        display(xo)
    end
end
#------------------------------------Fixed point method --------------------------------------

f = x -> (x.^3+4*x.^2).-10 # any function that a user wants
xo  = 0.5; # starting point or the fixed point
g =x-> ((x-2*x.^2).+3)^(1/4)
g = x -> ((2*x).-1+0im).^(1/3)
g =x-> 1/2*((x.^3).+1)

#------------P1-----------------------
#testing convergence ->
using ForwardDiff
using ChainRules
using Calculus
g = x->cos(x)
der_g =derivative(g)
x = -2*pi:0.0001:2*pi
der_g1 = der_g.(x)
for i in 1:length(x)
if(abs(der_g1[i])>=1.0)

display("It is divergent")
display(der_g1[i])
break;
end
if(i==length(x))
println("It is convergent")
display(der_g1[i])
end

end


#-------------------------P2----------------------
g = x-> exp(-x)
FixedPointMethod(g,xo)

#--------------------Fixed point Method---------
function FixedPointMethod(g,xo)
for i in 1:300
x_n1 = g.(xo)
if(abs(x_n1-xo)< 10^-5)
    display(x_n1)
    display(i)
    break;
end
xo = x_n1

end
end
#-----------------------------------------------------------------------------------------
#-----------------------------------Newton method --------------------------------------
using Plots
plotly();

f = x->x^2-7*x+10;
h = 10^-5
f1 = x-> (f.(x+h)-f(x))/h #derivative
xo = 1
#function NewtonMethod(f,f1,xo)
for i in 1:300
x_n1 = xo-f.(xo)*(f1.(xo)).^-1
    if(abs(x_n1-xo) < 10^-8)
            display(x_n1)
            display(["iterations: ", i])
            break;
    end
    xo = x_n1
    display(x_n1)

#end
end
x = -2*pi:0.0001:2*pi
y = x
plot!(x,g,linewidth=3,color="blue")
plot!(x,y,linewidth=3,color="red")


# plotting of a function ------------------------------
using Plots
x = -10:0.0001:10
plot(x,f,linewidth=3,color="blue")


#-------------------------problem 1----------------
using SymEngine
using ForwardDiff
using ChainRules
a = symbols(:a)
funt = 1/2+1/4*a.^2-a*sin.(a)-1/2*cos.(2*a)
f =x-> 1/2+1/4*x.^2-x*sin.(x)-1/2*cos.(2*x)
f1 = derivative(f)
xo  = pi/2
NewtonMethod(f,f1,xo)
# f(x) = 0 -> x = -+ 1.895494...and 0
xo = 5*pi
xo = 10*pi
x = -3:0.0001:3
plot(x,f,linewidth=3,color="blue")

f = x->(x-1)+2*x^3



#------------------------problem 2-------------------------------
x = -2:0.001:2
f =x-> x^2
#we need to find the "smallest" distance between some value of (x,y) and (1,0)
plot(x,f,linewidth=3,color="red")

#the distance between some (x,y) and (1,0) ->
d =x-> ((x-1)^2+(x^2)^2)^(1/2)
der_d = derivative(d)
f = der_d
f2 = derivative(f)
NewtonMethod.(f,f2,1.1) #anything in a region between -1000 and 1 should work
plot!(x,d,linewidth=3,color="green")
plot!(x,f,linewidth=3,color="blue")
plot!(x,f2,linewidth=3,color="black")



#------------------------------------Problem 3-------------------------
#sum of two numbers -> a+b =20
#(a+sqrt(a))*(b+sqrt(b)) = 155.55

#a = 20-b
x = 0:0.001:20
f = a-> (a+ a^(1/2))*((20-a)*(20-a)^(1/2))-155.55 #doesn't recognise it->
f =a-> (a+sqrt(a))*(sqrt(20-a))-a^(3/2)-a^2+20*sqrt(a)+20a-155.55 # matlab

plot(x,f,linewidth=3,color="blue")
NewtonMethod(f,derivative(f),5)
NewtonMethod(f,derivative(f),15)




#-----------------------------------problem 4---------------------------
x = 0:0.001:1
f = x-> 10^6*exp.(x)+4.35*10^5*x.^-1*(exp(x)-1)-1.564*10^6
NewtonMethod(f,derivative(f),1)
plot(x,f,linewidth=3,color="blue")


#-----------------------------------Newton method Use this One--------------------------------------
using Calculus
using ForwardDiff
using ChainRules
f = x-> x^10
f1 = x-> 10*x^9
xo = 1
function NewtonMethod(f,f1,xo)
for i in 1:200
x_n1 = xo-f.(xo)*(f1.(xo)).^-1
    if(abs(x_n1-xo) < 10^-8)
            display(x_n1)
            display(["iterations: ", i])
            return x_n1
            break;
    end
    xo = x_n1
    display(x_n1)

end
end
x = -2*pi:0.0001:2*pi
y = x
plot!(x,g,linewidth=3,color="blue")
plot!(x,y,linewidth=3,color="red")


# plotting of a function ------------------------------
using Plots
x = -10:0.0001:10
plot(x,f,linewidth=3,color="blue")





#-------------------------------------Secent Method -----------------------------------------------

f = x ->x^2-x-1;
x0 = -2;
x2= 0;

function SecantMethod(f,x1,xo)
for i in 1:2
display(x1)
x2 = x1-(f.(x1)*(x1-xo))*(f.(x1)-f.(xo)).^-1
xo = x1;
    if(abs(x2-x1)<10^-8)
        display(x2)
        display(i)
        break;
    end
x1 = x2;
end
end

#---------------------------Regula Falsi Method-----------------------
f = x ->x^2-x-1;
x0 = -2;
x2= 0;


function RegulaFalsi(f,x0,x2)

for i in 1:300

x1 = x0 - (f.(x0)*(x0-x2))*(f.(x0)-f.(x2)).^-1

    if(abs(x1-x0)<10^-8 || abs(x1-x2)<10^-8)
        display(x1)
        display(["number of iteratiosn : " ,i])
        break;

    end


    if(sign(f.(x1))*sign(f.(x0)) > 0 && sign(f.(x0)) < 0 && sign(f.(x2))>=0 )
        x0 = x1;
        display(x1)

    end

#                    Second case scenario
        if (sign(f.(x1))*sign(f.(x0)) > 0 && sign(f.(x0)) > 0 && sign(f.(x2))<=0)
        x0 = x1;
        display(x1)
    end

#                    Third case scenario
    if(sign(f.(x1))*sign(f.(x2)) > 0 && sign(f.(x2)) < 0 && sign(f.(x0)>=0))
        x2 = x1;
        display(x1)

    end
    #                Fourth case scenario
    if(sign(f.(x1))*sign(f.(x2)) > 0 && sign(f.(x2)) > 0 && sign(f.(x0)<=0))
        x2 = x1;
        display(x1)
    end
                #    Fifth case scenario
    if(f.(x1)==0)
    display(x1)
    end





end

end





#--------------------Fixed point Method---------
function FixedPointMethod(g,xo)
for i in 1:300
x_n1 = g.(xo)
if(abs(x_n1-xo)< 10^-5)
    display(x_n1)
    display(["Number of iterations :  ",i])
    break;
end
xo = x_n1

end
end

#---------------------Problem 1---------------------------------------------

f = x-> x^2-6
NewtonMethod(f,derivative(f),1)
SecantMethod(f,-1.5,2.5)
RegulaFalsi(f,-1.5,2.5)

f = x-> x^3-2*x^2-5
NewtonMethod(f,derivative(f),1)
SecantMethod(f,1,4)
RegulaFalsi(f,2,4)

f = x-> x-cos(x)
NewtonMethod(f,derivative(f),pi)
SecantMethod(f,1,4)
RegulaFalsi(f,1,4)

f = x-> x^3+3*x^2-1
NewtonMethod(f,derivative(f),2)
SecantMethod(f,1,4)
RegulaFalsi(f,0,4)

f = x-> x-0.8-0.2*sin(x)
NewtonMethod(f,derivative(f),2)
SecantMethod(f,1,4)
RegulaFalsi(f,0,4)

f = x-> 230*x^4+18*x^3+9*x^2-221*x-9
NewtonMethod(f,derivative(f),2)
SecantMethod(f,1,4)
RegulaFalsi(f,0.5,1)

#----------------------------------------------------------------------------
x = -2:0.001:2
y  = x
y1 = cos.(x)
plot!(x,y1,linewidth=3)
plot!(x,y,linewidth=3,color="red")


#applied Newton, Secant and RegulaFalsi

f = x-> (2*(x-1)+4*x^3)*(2*((x-1)^2+x^4)^(1/2)).^-1
NewtonMethod(f,derivative(f),0.5)

d = ((x-1)^2+x^4)^(1/2)

#-----------------------------------------------
f = x-> ((2x-4)-(2*(1/x-1)/(x^2)))/(2*(x-2)^2+(1/x-1)^2)^(1/2)
NewtonMethod(f,derivative(f),15)

f = b->((20-b)+(20-b)^(1/2))*(b+(b)^(1/2))-155.55

plot(0:0.0001:15,f,linewidth=3)

#---------------------------------------------------------------
f = l-> 10^6*exp.(2*l)+4.35*10^5/(l)*(exp.(l)-1)

plot(-5:0.00001:2,f,linewidth=3)

NewtonMethod(f,derivative(f),0.5)

f = i-> -1.35*10^5*i+1000*(1-(1+i)^(-360))
using Plots
plot(0:000001:1,f,linewidth=3)
SecantMethod(f,0.1,0.2)
BisectionMethod(f,0.0067,0.1)

RegulaFalsi(f,0,15)

f = x->x^2-3*x-1
NewtonMethod(f,derivative(f),0.3)

f = x->3^(-x)
SecantMethod(f,-1,0.5)

FixedPointMethod(f,-200)
f = x-> ((x+3)/(x^2+2))^(1/2)
f = x-> (2*x^2+1).^(1/3)

plot(2.1:0.0001:4,f,linewidth=3)

f =x-> pi+1/2*sin(x/2)

f = x->2+sin.(x)

FixedPointMethod(f,10000)

f = x->log(3*x^2)
FixedPointMethod(f,0.5)

f = x->cos.(x)
#--------------------------List problems Fixed point-------------------------#

f = x-> cos(x)
FixedPointMethod(f,0)

f  =x-> exp(-x)
FixedPointMethod(f,0)

f= x-> log(1+x)

f = x-> 1/2+1/16*(-1/2*log(1/2)-x*log(x))
FixedPointMethod(f,1/2)


f = x-> 1+exp(-x)





#----------------------------Problems ----------------------#
BisectionMethod(f,-1,2)
f =x-> 2*x+3*cos.(x)-exp.(x)
f = x->x*cos.(x) - 2x.^2 + 3*x - 1
f = x-> x+1-2*sin(pi*x)

f = x->x.^2-6
f1 = derivative(f);
xo = 1

x = -1:0.0001:1
using Plots
plot(x,f,linewidth=3)
plotly()

NewtonMethod(f,f1,-10)
SecantMethod(f,x1,xo)
FixedPointMethod(g,xo)

f =x-> x^3-25

f = x-> exp.(x)-2-cos.(exp.(x)-2)

f = h-> -12.4+10*(0.5*pi-asin.(h)-h*(1-h.^2).^(1/2))

plot(x,f,linewidth=3)

f = w->-1.7-32.17*(2*w.^2).^-1*((exp.(w)-exp.(-w))/2-sin.(w))

plot(x,f,linewidth=3)

f = x-> x.^2-6
xo = 3;
x1 = 2;
SecantMethod(f,0,0.48)
x = 0:0.0001:3
RegulaFalsi(f,0,0.48)

f = x->x-0.8-0.2*sin.(x)
xo = 3

f = x-> log(x-1)+cos.(x-1)
xo = 2
f = x-> (x-2).^2-log.(x)

f = x-> sin.(x)-exp(-x)
f = x-> (sin.(pi*x))/cos.(pi*x)-6

f = x-> 4*x.^2-exp.(x)-exp(-x)
f = x-> 2*(x-1)+4*x^3
NewtonMethod(f,derivative(f),0.5)
f1 = x->(x-2)^2+(1/x-1)^2
f = derivative(f)
f =i->-A*i+ P*(1-(1+i)^(-n))
A = 135000;
P = 1000;
n = 360;
x= 0:0.0001:0.9
f = x-> (1+2*x^2)^(1/3)
FixedPointMethod(f,1)
f = x-> (2*x-1)^(1/3)
f = x-> pi+1/2*sin(x/2)
f = x-> 2^(-x)
f = x-> 1/pi*asin(-x/2)+2
f = x-> -2*sin.(pi*x)
f = x->exp.(6*x)+3*log(2)*log(2)*exp.(2*x)-(log.(8))*exp(4*x)-(log(2)).^3
NewtonMethod(f,derivative(f),-1)
plot(-1:0.00001:1,f,linewidth=3)

R = 225;
C = 0.6*10^-6
L =0.5
Z = 100
f = w-> -1/Z + (R.^(-2)+(w*C- (w*L).^-1).^2)^(1/2)
plot(0:10:1000,f,linewidth=3)
BisectionMethod(f,10,250)
NewtonMethod(f,derivative(f),100)
SecantMethod(f,1,2)
RegulaFalsi(f,300,100)
y = 1
y0 = 0.8
g = 9.8
x = 90
v = 30
f = theta->sind.(theta)*cosd.(theta)*x-g*x.^2*(2*v.^2).^-1+y0*(cosd.(theta))^2
plot(0:0.001:3,f,linewidth=3)
BisectionMethod(f,0.11,6)
BisectionMethod(f,44,56)
NewtonMethod(f,derivative(f),5)
R = 3;
V = 30;
f = h-> -V+pi*h.^2*((3*R-h)/3)
f = h-> -V+(3*R*pi*h.^2)/3-(pi*h.^3)/3
f = x-> tan.(x)
FixedPointMethod(f,5)
f = x-> -1/x+(1/tan.(x))+x
f =x-> atan.(x)+pi
so = 300
m = 0.25
g = 32.17
k = 0.1
f = t->so-t*((m*g)/(k))+((m.^2*g)/(k.^2))*(1-exp.((-k.*t)/m))
FixedPointMethod(f,0)
NewtonMethod(f,derivative(f),5)
plot(0:0.001:10,f,linewidth=3)
f = x->x.^4-7*x.^3+18*x^2-20*x+8
NewtonMethod(f,derivative(f),-100)
f = x->14*x*exp.(x-2)-12*exp.(x-2)-7*x.^3+20*x^2-26*x+12
plot(0:0.0001:3,f,linewidth=3)
f = x-> exp.(x)+sin.(x)-4
f = x-> (x^3+1)/2
SecantMethod(f,1,2)
f = x-> (x^3+1)/2
f = x-> sin(x)+x
f1 = derivative(f)
g = derivative(f)

FixedPointMethod(f,2)

f = x->det([1 2 3 x;4 5 x 6;7 x 8 9;x 10 11 12])-1000
plot(-20:0.01:20,f,linewidth=3)
f = x-> x.^4-202*x^2+1404*x-3475
plot(-20:0.01:20,f,linewidth=3,color="red")
x = 1:1:10
Dump.(x)
Dump(1:1:10)





#-------------------------------------------------
f = x->x.^4-7x.^3+18*x.^2-20*x+8
g = x->8*x^4-12*x^3+6x^2-x
NewtonMethod(f,derivative(f),1)
NewtonMethod(g,derivative(g),3)
f = x-> 14*x*exp.(x-2)-12*exp(x-2)-7*x.^3+20*x.^2-26*x+12
f1 = derivative(f)
f2 = derivative(f1)
f3 = derivative(f2)



#--------------------Prepartion--------------------------------#
f =x-> x-tan.(x)
BisectionMethod(f,-2,0.78)
f = x->exp.(x)-2-cos.(exp.(x)-2)

x = -3:0.0001:3
plot(x,f,linewidth=3)

f = x-> (x+2).*(x+1).^2*x.*(x-1).^3*(x-2)

f = h-> -12.4+10*(0.5*3.14-asin.(h)-h*(1-h.^2)^(1/2))

f = w->-1.7*-32.17*((exp.(w)-exp.(-w))/2-sin.(w))

plot(x,f,linewidht=3)

f = x->(3+3*x^2)^(1/4)

FixedPointMethod(f,1)

NewtonMethod(f,derivative(f),5)
f = x->cos.(x)+x
FixedPointMethod(f,2)
x = -20:0.001:20
f = x-> det([1 2 3 x;4 5 x 6;7 x 8 9;x 10 11 12])-1000
plot(x,f,linewidth=3)
BisectionMethod(f,-20,-16)


f = x->(1+x).^(1/3)
FixedPointMethod(f,1)

f =x->pi+0.5*sin(x/2)
FixedPointMethod(f,1)


a = [1 -1 1;1 1 1;1 2 4]
f = [-6;0;6]
a1 = inv(a)*f

#--------------------------------------------------------
