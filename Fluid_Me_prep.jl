#--------------------------Prep--------------------------#
Re = [4000 10^4 10^5 10^6 10^7];
Le_over_d = 1.6*Re.^(1/4)

# fully rough flow


epsi_over_d = [0.00001 0.0001 0.001 0.01 0.05]
f = (-2*log10.((epsi_over_d)/3.7)).^-2


#---------------------------------Solving equation using NM---------#
epis_over_d = 0.0002
epis_over_d = 0;
Red = 1.29*10^5;
f = 0.01912
func = f-> (-1/f^(1/2))-2*log10(epis_over_d/3.7+(2.51)/(Red*f^(1/2)))
func1 = epis_over_d-> (-1/f^(1/2))-2*log10(epis_over_d/3.7+(2.51)/(Red*f^(1/2)))

h = 10^-5
xo = 0.01 #starting point - later on it turns out that 2 should have been taken as a starting point
function NewtonMethod(f,xo,h)
        f1 = x-> (f.(x+h)-f(x))/h; # approximation of derivative
        for i in 1:300
                x_n1 = xo-f.(xo)*(f1.(xo)).^-1

                    if(abs(x_n1-xo) < 10^-5)
                            display(x_n1)
                            display(["iterations: ", i])
                            break;
                    end
                    xo = x_n1
                    display(x_n1)

        end
end
#------------------------------------------------------#
xo = 0.000005
NewtonMethod(func1,xo,h)
NewtonMethod(func,xo,h);
h = 10^-5
xo = 0.1;
g = 9.81
d = 0.25*10^-3
L = 1.5*10^-2
h1 = 0.5
rho = 900
mju =0.002
func22 = v-> v^2/(2*g)+(-h1+(32*mju*L*v)/(rho*g*d^2))
NewtonMethod(func3,0.5,h);
func3 = v-> v^2+1800*v+0.1048
Q1 = 1.5*10^-7
A1 = d^2*pi/4
v = Q1/A1
hf = 32*(mju*v*L)/(rho*g*d^2)
fun5 = p-> p/(rho*g)+v^2/(2*g)+(hf)

NewtonMethod(fun5,5000,h);

A2 = (10^-2)^2*pi/4
F = -51138.69470513205*A2


y1 = v-> v^3+231-80*v
xo = 20
NewtonMethod(y1,xo,h)

f =0.013533252907878214
y2 = v-> 0.0509*v^2+3.058*v^2*f-4.2

 NewtonMethod(y2,2,h)
v1 = 6.122084732396374
Red = 10^3*v1*0.07/0.001
