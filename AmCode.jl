
using Plots
plotly();
#-------------------------------------------------Problem 1 ---------------------------------------#


#---------------------initial cond -----------------#
h_q = 6.5;
g = 32.2;
vo = 50;
xo = 60;
h_r = 7;

# x = teta
f = x-> xo*tan.(x)*(cos.(x))^2-1/2*((xo^2*g)/(vo^2))-(h_r-h_q)*(cos.(x))^2
a = 0;
b = 1;

function BisectionMethod(f,a,b)

    epsi = 10^-10

        n = (log(abs(a-b)/epsi)).*log(2)^-1
        n = round(n,digits=0)

    #               finding the roots
    for i in 1:n

    bis = a+((b-a)/2);
    f_bis = f.(bis);
    f_a = f.(a);
    f_b = f.(b);

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
 BisectionMethod(f,a,b)
# angle is 0.45384738000575453 -> in degrees
deg = (0.45384738000575453/3.14)*180
x1 = -1:0.001:1
plot(x1,f.(x1),linewidth="3")






#------------------------------------------Problem 2------------------------------------#

#---------------------initial cond -----------------#
Q = 9.4*10^-6
q = 2.4*10^-5
R = 0.1
F =0.3
e_0 = 0.885*10^-12

#-------------------------------------------------Problem 1 ---------------------------------------#

f = z->-0.3+((Q*q*z)/(2*e_0))*(1-z./((z.^2+R^2)^(1/2)))
a = 0;
b = 1;
BisectionMethod(f,a,b)
# result is 0.002411878085695207 -> epsi = 10^-10
# does the result make sense? with such a miniature force 2.4mm distance is reasonable
x1 = -1:0.00001:1;
plot(x1,f.(x1),linewidth="3")
# it crosses the x axis at the given point -> QED

#--------------------Problem 3------------------------------#

#---------------------Initial conditions----------------------#

wo = 20*10^3
I = 52.9*10^-6
E = 70*10^9
L = 4
C1 = (wo)/(360*L*E*I); #to write it in a compact way

#------------------------------------------------------------#
f = x-> C1*(7*L^4-10*L^2*x^2+3*x^4)+C1*x*(-20*L^2*x+12*x^3)
h = 10^-5
f2 =x-> (1/(x^(1/2)))+2*log10((0.000125/3.7)+(2.51/(433000*x^(1/2))))
xo = 0.5 #starting point - later on it turns out that 2 should have been taken as a starting point


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

NewtonMethod(f,xo,h);
# result 5.261629775415293 (there are three other roots)
x1 = -1:0.00001:1;
plot(x1,f.(x1),linewidth="3")

# when we find the root x for the max deflection -> put it back in the original eq.
f1 = x->C1*x*(7*L^4-10*L^2*x^2+3*x^4)
maxDef = abs(f1(5.261629775415293));

# Secant method
f = x-> C1*(7*L^4-10*L^2*x^2+3*x^4)+C1*x*(-20*L^2*x+12*x^3)
x0 = 1.5;
x1= 2.5;

function SecantMethod(f,xo,x1)
    for i in 1:100
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

SecantMethod(f,x0,x1);
#result ->  2.0773184894369123
maxDef_new = f1(2.0773184894369123) #0.0090179

#why are results different? The nonlinear function is of the 4th order and has 4 roots -> thus
#our first assumption that x = 4 does not produce the maximum deflection, it turns out that if we
#were to start from 2, we would have reached a better approximation of the maximum deflection.


#----------------------------------Problem 4----------------------------#
# x1 = Tc , x2 = J_c, x3 = T_h, x4 = J_h

# Jacobian is defined as 
#                _                                                                         _
#               |   dy1/dx1      dy1/dx2    dy1/dx3    . . .                   dy1/dx_n     |
#               |   dy2/dx1      dy2/dx2    dy2/dx3                                         |
#               |   dy3/dx1      dy3/dx2    dy4/dx3 .                                       |
#               |   dy4/dx1      dy4/dx2    dy4/dx3     .                                   |
#       J_f =   |   dy5/dx1         .          .            .                               |
#               |   dy6/dx1         .          .               .                            |
#               |   dy7/dx1         .          .                    .                       |
#               |      .            .          .                        .                   |
#               |      .            .          .                            .               |
#               |   dy_n/dx1     dy_n/dx2   dy_n/dx3    . . .                  dy_n/dx_n    |
#                --                                                                       --
#
#
#



df_p(x1,x2,x3,x4)= [5.67*10^-8*x1^4+17.41*x1-x2;x2-0.71*x4+7.46*x1;5.67*10^-8*x3^4+1.865*x3-x4;x4-0.71*x2+7.46*x3];

df_dx1(x1,x2,x3,x4) = (df_p.(x1+h,x2,x3,x4)-df_p(x1,x2,x3,x4))/h
df_dx2(x1,x2,x3,x4) = (df_p.(x1,x2+h,x3,x4)-df_p(x1,x2,x3,x4))/h
df_dx3(x1,x2,x3,x4) = (df_p.(x1,x2,x3+h,x4)-df_p(x1,x2,x3,x4))/h
df_dx4(x1,x2,x3,x4) = (df_p.(x1,x2,x3,x4+h)-df_p(x1,x2,x3,x4))/h

h = 10^-5;
J = [df_dx1;df_dx2;df_dx3;df_dx4];
J2(x1,x2,x3,x4) = [J[1](x1,x2,x3,x4) J[2](x1,x2,x3,x4) J[3](x1,x2,x3,x4) J[4](x1,x2,x3,x4)]


x_i =[298*10^3;3000;298*10^3;5000];
tol = 10^-6
counter = 0;
y_initial = [5188.18;2352.71;2250;11093];


for i in 1:100

    x_t = ((J2.(x_i[1],x_i[2],x_i[3],x_i[4]))\(df_p.(x_i[1],x_i[2],x_i[3],x_i[4])-y_initial))
    x_next = x_i-x_t
    display(x_next);
    for j in 1:length(x_i)
        if(abs(x_next[j]-x_i[j])<=tol)
            display(["Iterations : ",i])
            display(x_next);
            counter = counter+1;
        end
        if(counter!=length(x_i) && j==length(x_i))
            counter = 0;
        end
    end
    x_i = x_next;
end

df_p(x1)= [5.67*10^-8*x1[1]^4+17.41*x1[1]-x1[2];x1[2]-0.71*x1[4]+7.46*x1[1];5.67*10^-8*x1[3]^4+1.865*x1[3]-x1[4];x1[4]-0.71*x1[2]+7.46*x1[3]];


function Jacobain(df_p,x1,numberOfVariables,h) # (matrix, inputs, number of Variables (x1,x2 ... xn))
    
    j = zeros(1,numberOfVariables)
    j_standard = zeros(1,numberOfVariables)

    for i1 in 1:numberOfVariables 

        j_standard[i1] = x1[i1];  
        j[i1] = x1[i1]; 

    end 
    J = zeros(1,numberOfVariables)' # default value 

    for i in 1:numberOfVariables  
    
    j[i] = j[i]+h
    J = [J (df_p(j)-df_p(j_standard))/h]
    display(J)
    j[i] = j_standard[i]; 

    end 
    
    J = J[:,2:size(J)[2]] #nullfies the input zero from the storage 

    return J

end 


J11 = Jacobain(df_p,[1 2 3 4],4, 10^-5)

(df_p([1+10^-5 2 3 4])-df_p([1 2 3 4]))/10^-5


# results ->
# Tc = 481.0272546993752
#Jc = 6222.225082461532
#Jh = 10504.19493312517
#Th = 671.1239779386751




#  problem 25 

f(t) = log(t) 
h = 10^-5
f_der(t) = (f(t+h)-f(t))/h 

c_prime_norm(t) = (1+t^(-2))^(1/2)

T1(t) = (1/c_prime_norm(t))
T2(t) = (t^-1*(c_prime_norm(t))^-1)

f_der1(t) = (T1(t+h)-T1(t))/h 
f_der2(t) = (T2(t+h)-T2(t))/h 




k(t) = (((f_der1.(t)).^2+(f_der2.(t)).^2)^(1/2))/c_prime_norm.(t)



f_der3(t) = (k.(t+h)-k.(t))/h 



NewtonMethod(f_der3,0.5,10^-5) # gives 0.7071

k(0.7071)  # 0.3848983650052994


