#----------------------Gaus Seidel Method -> one function-----------##

#using Plots;
#plotly();
using SymPy;

a= symbols("a");
b = symbols("b");

f = x-> x.^3-6*x^2+9*x-4;
i = 9; # number in front of the x ->
f1 =f(a);
x_new =-i.^-1*(f(b)-i*b);

x_initial = 3.5;
tol = 10^-6;
for j in 1:100

    x_new1 = x_new.(x_initial);
    display(x_new1);
    if(abs((x_initial-x_new1))<=tol)
        display(["x :",x_new1])
        display(["Iterations : ",j]);
        break;
    end
    x_initial = x_new1;


end

#-----------------Gaus Seidel Method for linear systems---------------

# in julia imaginary number => im = sqrt(-1)

A = [10-4im -2-10im -4+5im;-2 6+3im -2;-4 -2 10+20im];
L = [10-4im -2-10im -4+5im;-2 6+3im -2;-4 -2 10+20im];
U = [10-4im -2-10im -4+5im;-2 6+3im -2;-4 -2 10+20im];
b = [-2;3;-1];

for i in 1 :length(A[1,:])
    for m in 1:length(A[1,:])
        if(m>i)
            L[i,m] = 0;
        end
    end
end


for i in 1 :length(A[1,:])
    for m in 1:length(A[1,:])
        if(m<i || m==i)
            U[i,m] = 0;
        end
    end

end

x_initial = [7;3;10];
tol = 10^-6;
counter = 0;

T = -L\U;
R = L\b;

for j in 1:150

    x_next = T*x_initial+R
    display(x_next)
    for i in 1:length(x_initial)
    if(abs((x_next[i]-x_initial[i]))<=tol)
        display(["Iterations : ",j])
        display(x_next);
        counter = counter+1;
    end
    if(i==length(x_initial) && counter!=length(x_initial))
        counter = 0;
    end
    end
    x_initial = x_next;
    if(counter == length(x_initial))
        display(["Iterations : ",j])
        display(x_next);
        break;
    end

end


 #---------------------Gauss Siedel method -> system without using matrix---

V_i1 = 1.05 # initial point
V_i2 = 1
V_i3 = 1;
tol = 10^-6
Y = [20-50im -(10-20im) -(10-30im);-10+20im 26-52im  -16+32im;-10+30im -16+32im 26-62im]
P2 =-2.566
Q2 = -1.102im
P3 =-1.386
Q3 =-0.452im
r = 2;
k = 1;
ynew =0;

for i in 1:20

    ynew = 0;
    for m = 1:length(Y[:,1])
        if(r!=k)
            ynew = ynew+Y[r,k]
        end
        k = k+1;
    end
    r = 2
    k = 1;
    V_next2 = (-1/ynew)*((P2-Q2)/conj(V_i2)-(Y[2,1]*V_i1+Y[2,3]*V_i3))
    r = r+1;
    ynew = 0
    for m = 1:length(Y[:,1])
        if(r!=k)
            ynew = ynew+Y[r,k]
        end
        k = k+1;
    end
    V_next3 = (-1/ynew)*((P3-Q3)/conj(V_i3)-(Y[3,1]*V_i1+Y[3,2]*V_next2))
    if(abs(V_next2-V_i2)<=tol)
        display([V_next2 V_next3])
        display(["Iterations : ",i])
        break;
    end
    V_i2 = V_next2;
    V_i3 = V_next3;
    display([V_next2 V_next3])
    r = 2
    k = 1;
end
P1_Q1 = -conj(V_i1)*(V_i1*(Y[1,2]+Y[1,3])-(Y[1,2]*V_i2+Y[1,3]*V_i3))
P1 = real(P1_Q1)
Q1 = imag(P1_Q1)
#-------------------Newton - Raphson's method -> one function---------------------#



f = x->x.^4+3*x.^3-15*x.^2-19*x+30
x_i =4;
h = 10^-6
f3 =x->(f.(x+h)-f.(x))*h.^-1
tol = 10^-6
y_initial = 7;
for i in 1:150
    x_t = ((f3.(x_i)))\(f.(x_i)-y_initial)
    x_next = x_i-x_t
    display(x_next);
    if(abs(x_next-x_i)<=tol)
        display(["Iterations : ",i])
        display(x_next);
        break;
    end
    x_i = x_next;
end


#---------------------Newton-Raphson's with Jacobian -> system (2x2)
using ForwardDiff;
df_p(x1,x2)= [x1^2+x2^2;exp.(x1)+x2];
df_dx(x1,x2) = ForwardDiff.derivative(x1->df_p(x1,x2),x1)
df_dy(x1,x2) = ForwardDiff.derivative(x2->df_p(x1,x2),x2)

J = [df_dx;df_dy];
J2(x1,x2) = [J[1](x1,x2) J[2](x1,x2)]

# J[1](1,1) -> calling it with proper ref ->


x_i =[1;1];
tol = 10^-6
counter = 0;
y_initial = [4;1];

for i in 1:100
    x_next = x_i - ((J2.(x_i[1],x_i[2])).^-1*(df_p.(x_i[1],x_i[2])-y_initial))
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
    if(counter==length(x_i))
        display(["Iterations : ",i])
        display(x_next);
        break;
    end
end

# partial derivatives

g(x,y) = [2*y^2+2*x^2+x*y,x^2+y^2]
dg_dx(x,y) = ForwardDiff.derivative(x->g(x,y),x)
dg_dy(x,y) = ForwardDiff.derivative(y->g(x,y),y)

#---------------------Newton-Raphson's with Jacobian -> system (3x3)
using ForwardDiff;
df_p(x1,x2,x3)= [x1^2-x2^2+x3^2;x1*x2+x2^2-3*x3;x1-x1*x3+x2*x3];
df_dx(x1,x2,x3) = ForwardDiff.derivative(x1->df_p(x1,x2,x3),x1)
df_dy(x1,x2,x3) = ForwardDiff.derivative(x2->df_p(x1,x2,x3),x2)
df_dz(x1,x2,x3) = ForwardDiff.derivative(x3->df_p(x1,x2,x3),x3)

J = [df_dx;df_dy;df_dz];
J2(x1,x2,x3) = [J[1](x1,x2,x3) J[2](x1,x2,x3) J[3](x1,x2,x3)]

# J[1](1,1) -> calling it with proper ref ->


x_i =[1;1;1];
tol = 10^-6
counter = 0;
y_initial = [11;3;6];

for i in 1:100
    # avoid this implementation use A\b not INVERSE
    x_next = x_i - ((J2.(x_i[1],x_i[2],x_i[3])).^-1*(df_p.(x_i[1],x_i[2],x_i[3])-y_initial))
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
    if(counter==length(x_i))
        display(["Iterations : ",i])
        display(x_next);
        break;
    end
end


#----------------------Numerical NR with partial derivaties--------------------#
df_p(x1,x2,x3)= [x1^2-x2^2+x3^2;x1*x2+x2^2-3*x3;x1-x1*x3+x2*x3];
h = 10^-6;
df_dx1(x1,x2,x3) = (df_p.(x1+h,x2,x3)-df_p(x1,x2,x3))/h
df_dx2(x1,x2,x3) = (df_p.(x1,x2+h,x3)-df_p(x1,x2,x3))/h
df_dx3(x1,x2,x3) = (df_p.(x1,x2,x3+h)-df_p(x1,x2,x3))/h

J = [df_dx1;df_dx2;df_dx3];
J2(x1,x2,x3) = [J[1](x1,x2,x3) J[2](x1,x2,x3) J[3](x1,x2,x3)]

# J[1](1,1) -> calling it with proper ref ->


x_i =[1;1;1];
tol = 10^-6
counter = 0;
y_initial = [11;3;6];

for i in 1:100
    x_t = ((J2.(x_i[1],x_i[2],x_i[3]))\(df_p.(x_i[1],x_i[2],x_i[3])-y_initial))
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
    if(counter==length(x_i))
        display(["Iterations : ",i])
        display(x_next);
        break;
    end
end



#-------------------------------NR with Numerical two variables
df_p(x1,x2)= [x2-3*x1+1.9;x2+x1.^2-1.8];
h = 10^-6;
df_dx1(x1,x2) = (df_p.(x1+h,x2)-df_p(x1,x2))/h
df_dx2(x1,x2) = (df_p.(x1,x2+h)-df_p(x1,x2))/h


J = [df_dx1;df_dx2];
J2(x1,x2) = [J[1](x1,x2) J[2](x1,x2)]

# J[1](1,1) -> calling it with proper ref ->


x_i =[1;1];
tol = 10^-6
counter = 0;
y_initial = [0;0];

for i in 1:100
    x_t = ((J2.(x_i[1],x_i[2]))\(df_p.(x_i[1],x_i[2])-y_initial))
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
    if(counter==length(x_i))
        display(["Iterations : ",i])
        display(x_next);
        break;
    end
end


#----------------------------Problem 28-------------------------------
Y = [5-10im -2+4im -3+6im;-2+4im 2-4im 0; -3+6im 0 3-6im]
V = [1;]
k = 2;

for i in 1:length(Y[:,1])

    Pi = V[k] *Y[k,i]V[i]

end

##----------------------------------------------

V_i1 = 1.0 # initial point
V_i2 = 1
tol = 10^-6
Y = [(4-8im) -(4-8im);-(4-8im) (4-8im)  ]
P1 =-0.7
Q1 = -0.4im
P2 =-0.7
Q2 = -0.4im
r = 1;
k = 1;
ynew =0;

for i in 1:1

    ynew = 0;
    for m = 1:length(Y[:,1])
        if(r==k)
            ynew = ynew+Y[r,k]
        end
        k = k+1;
    end
    k = 1;
    V_next1 = (-1/ynew)*((P1-Q1)/conj(V_i1)-(Y[1,2]*V_i2))
    r = r+1;
    ynew = 0
    for m = 1:length(Y[:,1])
        if(r==k)
            ynew = ynew+Y[r,k]
        end
        k = k+1;
    end
    V_next2 = (-1/ynew)*((P2-Q1)/conj(V_i2)-(Y[2,1]*V_next1))
    if(abs(V_next2-V_i2)<=tol)
        display([V_next1 V_next2])
        display(["Iterations : ",i])
        break;
    end
    V_i2 = V_next2;
    V_i1 = V_next1
    display([V_next1 V_next2])
    r = 1;
    k = 1;
end
P1_Q1 = -conj(V_i1)*(V_i1*(Y[1,2]+Y[1,3])-(Y[1,2]*V_i2+Y[1,3]*V_i3))
P1 = real(P1_Q1)
Q1 = imag(P1_Q1)




# ------------Gaus Siedel method with PS------------#
V_i1 = 0.866+0.5im;
V_i2 = 1.0;

V = [V_i1 V_i2];

P1 = -0.7;
Q1 = -0.4im;
P2 =-0.7;
Q2 = -0.4im;
tol = 10^-6
Y = [(4-8im) -(4-8im);-(4-8im) (4-8im)  ]
ynew = 0;

R = 2;
r = 1;
k = 1;
Ykn = 0;

for j in 1:length(Y[:,1])

    if(r==k)
        ynew = ynew + Y[r,k];
    end
    k = k+1;
end


for i in 1:30

        k = 1;
        Ykn = 0;
        for m in 1:length(Y[:,1])
            if(R!=m)
            Ykn = Ykn+Y[R,m]*V[m];
        end
        end

    V_next2 = 1/ynew*((P2-Q2)/(conj(V[R]))-(Ykn));
    if(abs(V_next2-V_i2)<=tol)
        display(["Iterations : ",i])
        display(V_next2);
        break;
    end
    display(V_next2);
    V_i2= (V_next2);
    V = [V_i1 V_i2]
end


#---------------------------Problem 34---------------------------#
V_i1 = 1.0;
V_i2 = 1.0; #guess
V_i3 =1 #guesses
V = [V_i1 V_i2 V_i3];

P2 =1
P3 = 1.5
Q2 = 0.5im
Q3 = 0.75im
tol = 10^-8
Y = [-10im 5im 5im;5im -10im 5im;5im 2im -10im ]
ynew = 0;

R = 2;
r = 2;
k = 1;
Ykn = 0;

for j in 1:length(Y[:,1])

    if(r==k)
        ynew = ynew + Y[r,k];
    end
    k = k+1;
end

r = 3
k = 1;
ynew1 = 0;
for j in 1:length(Y[:,1])

    if(r==k)
        ynew1 = ynew1 + Y[r,k];
    end
    k = k+1;
end

    for i in 1:150
        R = 2;
        k = 1;
        Ykn = 0;
        for m in 1:length(Y[:,1])
            if(R!=m)
            Ykn = Ykn+Y[R,m]*V[m];
        end
        end


    V_next2 = 1/ynew*((P2-Q2)/(conj(V[R]))-(Ykn));
    R = R+1;
    Ykn = 0;
    V_i2= (V_next2);
    V = [V_i1 V_i2 V_i3];

    for m in 1:length(Y[:,1])
        if(R!=m)
        Ykn = Ykn+Y[R,m]*V[m];
    end
    end
    R = 2

    V_next3 = 1/ynew1*((P3-Q3)/(conj(V_i3))-(Ykn));
    if(abs(V_next2-V_i2)<=tol && abs(V_next3-V_i3)<=tol)
        display(["Iterations : ",i])
        display([V_next2,V_next3]);
        break;
    end
    display(["Iterations : ",i])
    display([V_next2,V_next3]);
    V_i3= (V_next3);
    V = [V_i1 V_i2 V_i3];
end



#------------------------------Problem 36----------------------#

V_i1 = 1.04;
V_i2 =1;
V_i3 = 1;
V_i4 = 1;
P2 =-0.5
Q2 =0.2im
P3 =1
Q3 =-0.5im
P4 =-0.3
Q4 =0.1im
tol = 10^-6
Y = [3-9im -2+6im -1+3im 0; -2+6im 3.666-11im -0.666+2im -1+3im;
    -1+3im -0.666+2im 3.666-11im -2+6im;0 -1+3im -2+6im 3-9im]
count = 0;
for i in 1:250


V_next2 = (1/Y[2,2])*((P2-Q2)/(conj(V_i2))-(Y[2,1]*V_i1+Y[2,3]*V_i3+Y[2,4]*V_i4));
display(["Iterations : ",i])
display([V_next2]);
if(abs(V_next2-V_i2)<=tol)
    display(["Iterations : ",i])
    display([V_next2]);
    count = count+1;

end

V_i2 = V_next2;
V_next3 = (1/Y[3,3])*((P3-Q3)/(conj(V_i3))-(Y[3,1]*V_i1+Y[3,2]*V_i2+Y[3,4]*V_i4));

if(abs(V_next3-V_i3)<=tol)
    display(["Iterations : ",i])
    display([V_next3]);
    count = count+1;
end

V_i3 = V_next3;
V_next4 = (1/Y[4,4])*((P4-Q4)/(conj(V_i4))-(Y[4,1]*V_i1+Y[4,2]*V_i2+Y[4,3]*V_i3));
if(abs((V_next4-V_i4))<=tol)
    display(["Iterations : ",i])
    display([V_next4]);
    count = count+1;
end
V_i4 = V_next4;
if(count ==3)
display("All of them converged")
break;
count = 0
end
end

#-------------------------------Newton Raphson method problem 38 --------------------------#

V_i1 = 1
V_i2 = 1
V_i3 = 1;
teta1 =0
teta2 = 0*pi/180;
teta3 =0*pi/180;
h = 10^-6;
Y= [-10im 5im 5im;5im -10im 5im;5im 5im -10im];

#P2 = V_i2*(abs(Y[2,1])*V_i1*cos(teta2-teta1-atan(imag(Y[2,1])/real(Y[2,1])))+abs(Y[2,2])*V_i2*cos(teta2-teta2-atan(imag(Y[2,2])/real(Y[2,2]))))
#2 = V_i2*(abs(Y[2,1])*V_i1*sin(teta2-teta1-atan(imag(Y[2,1])/real(Y[2,1])))+abs(Y[2,1])*V_i2*sin(-atan(imag(Y[2,2])/real(Y[2,2]))))


#dP2(V_i2,teta2) = V_i2*((abs(Y[2,1])*V_i1*cosd(teta2-teta1-atand.(imag(Y[2,1])/real(Y[2,1]))))+
#(abs(Y[2,2])*V_i2*cosd.(-atand(imag(Y[2,2])/real(Y[2,2]))))+(abs(Y[2,3])*V_i3*cosd.(teta2-teta3-atand(imag(Y[2,3])/real(Y[2,3])))))

#dP2(V_i2,teta2) = -V_i2*((abs(Y[2,1])*V_i1*sind(teta2-teta1-atand.(imag(Y[2,1])/real(Y[2,1]))))+
#(abs(Y[2,3])*V_i3*sind.(teta2-teta3-atand(imag(Y[2,3])/real(Y[2,3])))))


# Jacobian for the partials of P2

dP2(V_i2,V_i3,teta2,teta3) = V_i2*((abs(Y[2,1])*V_i1*cos(teta2-teta1-atan.(imag(Y[2,1])/real(Y[2,1]))))+
(abs(Y[2,2])*V_i2*cos.(-atan(imag(Y[2,2])/real(Y[2,2]))))+(abs(Y[2,3])*V_i3*cos.(teta2-teta3-atan(imag(Y[2,3])/real(Y[2,3])))))

display("------------DP2--------")

dP2_dteta22(V_i2,V_i3,teta2,teta3) =(dP2.(V_i2,V_i3,teta2+h,teta3)-dP2.(V_i2,V_i3,teta2,teta3))/h
dP2_dV22(V_i2,V_i3,teta2,teta3)=(dP2.(V_i2+h,V_i3,teta2,teta3)-dP2.(V_i2,V_i3,teta2,teta3))/h
dP22_dteta3(V_i2,V_i3,teta2,teta3)=(dP2.(V_i2,V_i3,teta2,teta3+h)-dP2.(V_i2,V_i3,teta2,teta3))/h
dP22_dV3(V_i2,V_i3,teta2,teta3)=(dP2.(V_i2,V_i3+h,teta2,teta3)-dP2.(V_i2,V_i3,teta2,teta3))/h

display(dP2_dteta22(V_i2,V_i3,teta2,teta3))
display(dP2_dV22(V_i2,V_i3,teta2,teta3))
display(dP22_dteta3(V_i2,V_i3,teta2,teta3))
display(dP22_dV3(V_i2,V_i3,teta2,teta3))

#Jacobian for the partials of Q2

dQ2(V_i2,V_i3,teta2,teta3) = V_i2*((abs(Y[2,1])*V_i1*sin(teta2-teta1-atan.(imag(Y[2,1])/real(Y[2,1]))))+
(abs(Y[2,2])*V_i2*sin.(-atan(imag(Y[2,2])/real(Y[2,2]))))+(abs(Y[2,3])*V_i3*sin.(teta2-teta3-atan(imag(Y[2,3])/real(Y[2,3])))))

dQ2_dteta22(V_i2,V_i3,teta2,teta3)=(dQ2.(V_i2,V_i3,teta2+h,teta3)-dQ2.(V_i2,V_i3,teta2,teta3))/h
dQ2_dV22(V_i2,V_i3,teta2,teta3) =(dQ2.(V_i2+h,V_i3,teta2,teta3)-dQ2.(V_i2,V_i3,teta2,teta3))/h
dQ3_dteta3(V_i2,V_i3,teta2,teta3)=(dQ2.(V_i2,V_i3,teta2,teta3+h)-dQ2.(V_i2,V_i3,teta2,teta3))/h
dQ3_dV3(V_i2,V_i3,teta2,teta3)=(dQ2.(V_i2,V_i3+h,teta2,teta3)-dQ2.(V_i2,V_i3,teta2,teta3))/h

display("------------DQ2--------")
display(dQ2_dteta22(V_i2,V_i3,teta2,teta3))
display(dQ2_dV22(V_i2,V_i3,teta2,teta3))
display(dQ3_dteta3(V_i2,V_i3,teta2,teta3))
display(dQ3_dV3(V_i2,V_i3,teta2,teta3))


#Jacobian for the P3

k = 3
r = 1;
dP3(V_i2,V_i3,teta2,teta3) = V_i3*((abs(Y[k,r])*V_i1*cos(teta3-teta1-atan.(imag(Y[k,r])/real(Y[k,r]))))+
(abs(Y[k,r+1])*V_i2*cos.(teta3-teta2-atan(imag(Y[k,r+1])/real(Y[k,r+1]))))+(abs(Y[k,r+2])*V_i3*cos.(teta3-teta3-atan(imag(Y[k,r+2])/real(Y[k,r+2])))))

display("------------DP3--------")

dP33_dteta22(V_i2,V_i3,teta2,teta3) =(dP3.(V_i2,V_i3,teta2+h,teta3)-dP3.(V_i2,V_i3,teta2,teta3))/h
dP33_dV22(V_i2,V_i3,teta2,teta3) =(dP3.(V_i2+h,V_i3,teta2,teta3)-dP3.(V_i2,V_i3,teta2,teta3))/h
dP33_dteta3(V_i2,V_i3,teta2,teta3)=(dP3.(V_i2,V_i3,teta2,teta3+h)-dP3.(V_i2,V_i3,teta2,teta3))/h
dP33_dV3(V_i2,V_i3,teta2,teta3)=(dP3.(V_i2,V_i3+h,teta2,teta3)-dP3.(V_i2,V_i3,teta2,teta3))/h

display(dP33_dteta22(V_i2,V_i3,teta2,teta3))
display(dP33_dteta3(V_i2,V_i3,teta2,teta3))
display(dP33_dV22(V_i2,V_i3,teta2,teta3))
display(dP33_dV3(V_i2,V_i3,teta2,teta3))

# Jacobian for the Q

dQ3(V_i2,V_i3,teta2,teta3) = V_i3*((abs(Y[k,r])*V_i1*sin(teta3-teta1-atan.(imag(Y[k,r])/real(Y[k,r]))))+
(abs(Y[k,r+1])*V_i2*sin.(teta3-teta2-atan(imag(Y[k,r+1])/real(Y[k,r+1]))))+(abs(Y[k,r+2])*V_i3*sin.(teta3-teta3-atan(imag(Y[k,r+2])/real(Y[k,r+2])))))


dQ33_dteta22(V_i2,V_i3,teta2,teta3) =(dQ3.(V_i2,V_i3,teta2+h,teta3)-dQ3.(V_i2,V_i3,teta2,teta3))/h
dQ33_dV22(V_i2,V_i3,teta2,teta3) =(dQ3.(V_i2+h,V_i3,teta2,teta3)-dQ3.(V_i2,V_i3,teta2,teta3))/h
dQ33_dteta3(V_i2,V_i3,teta2,teta3)=(dQ3.(V_i2,V_i3,teta2,teta3+h)-dQ3.(V_i2,V_i3,teta2,teta3))/h
dQ33_dV3(V_i2,V_i3,teta2,teta3)=(dQ3.(V_i2,V_i3+h,teta2,teta3)-dQ3.(V_i2,V_i3,teta2,teta3))/h
display("------------DQ3--------")

display(dQ33_dteta22(V_i2,V_i3,teta2,teta3))
display(dQ33_dteta3(V_i2,V_i3,teta2,teta3))
display(dQ33_dV22(V_i2,V_i3,teta2,teta3))
display(dQ33_dV3(V_i2,V_i3,teta2,teta3))


J5(V_i2,V_i3,teta2,teta3) = [dP2_dteta22(V_i2,V_i3,teta2,teta3) dP22_dteta3(V_i2,V_i3,teta2,teta3) dP2_dV22(V_i2,V_i3,teta2,teta3) dP22_dV3(V_i2,V_i3,teta2,teta3);
dP33_dteta22(V_i2,V_i3,teta2,teta3) dP33_dteta3(V_i2,V_i3,teta2,teta3) dP33_dV22(V_i2,V_i3,teta2,teta3) dP33_dV3(V_i2,V_i3,teta2,teta3) ;
dQ2_dteta22(V_i2,V_i3,teta2,teta3) dP33_dV22(V_i2,V_i3,teta2,teta3) dQ2_dV22(V_i2,V_i3,teta2,teta3)  dQ3_dV3(V_i2,V_i3,teta2,teta3);
dQ33_dteta22(V_i2,V_i3,teta2,teta3) dQ33_dteta3(V_i2,V_i3,teta2,teta3) dQ33_dV22(V_i2,V_i3,teta2,teta3)  dQ33_dV3(V_i2,V_i3,teta2,teta3)
]


y_initial1(V_i2,V_i3,teta2,teta3) = [dP2.(V_i2,V_i3,teta2,teta3);dP3.(V_i2,V_i3,teta2,teta3);dQ2.(V_i2,V_i3,teta2,teta3);dQ3.(V_i2,V_i3,teta2,teta3)]

P2 =-1
P3 = -1.5
Q2 = -0.5im
Q3 = -0.75im

y = [P2;P3;Q2;Q3]
delY = y_initial1(V_i2,V_i3,teta2,teta3)-y; # for my NR iterations


tol = 10^-2
counter = 0;
X_in = [teta2;teta3;V_i2;V_i3]
for i in 1:100

    X_t = ((J5.(X_in[3],X_in[4],X_in[1],X_in[2])\(delY)))
    X_next = X_in-X_t # book does it in an inverse way keep that in mind
    display(X_next[3])
    display("------------")
    display(X_next[2])
        if(abs(X_next[3]-X_in[3])<=tol)
            display(["Iterations : ",i])
            display(X_next[3]);
            break;
        end
    X_in =X_next;
    delY = y_initial1(V_i2,V_i3,teta2,teta3)-y; # for my NR iterations
end



a = dP2(X_in[3],X_in[4],X_in[1],X_in[2])





b = J5.(X_in[3],X_in[4],X_in[1],X_in[2])



#-------------------------------Newton Raphson method problem 43 --------------------------#

V_i1 = 1
V_i2 = 1
V_i3 = 1;
teta1 =0
teta2 = 0*pi/180;
teta3 =0*pi/180;
h = 10^-6;
Y= [-12.5im 10im 2.5im;10im -15im 5im;2.5im 5im -7.5im];

#P2 = V_i2*(abs(Y[2,1])*V_i1*cos(teta2-teta1-atan(imag(Y[2,1])/real(Y[2,1])))+abs(Y[2,2])*V_i2*cos(teta2-teta2-atan(imag(Y[2,2])/real(Y[2,2]))))
#2 = V_i2*(abs(Y[2,1])*V_i1*sin(teta2-teta1-atan(imag(Y[2,1])/real(Y[2,1])))+abs(Y[2,1])*V_i2*sin(-atan(imag(Y[2,2])/real(Y[2,2]))))


#dP2(V_i2,teta2) = V_i2*((abs(Y[2,1])*V_i1*cosd(teta2-teta1-atand.(imag(Y[2,1])/real(Y[2,1]))))+
#(abs(Y[2,2])*V_i2*cosd.(-atand(imag(Y[2,2])/real(Y[2,2]))))+(abs(Y[2,3])*V_i3*cosd.(teta2-teta3-atand(imag(Y[2,3])/real(Y[2,3])))))

#dP2(V_i2,teta2) = -V_i2*((abs(Y[2,1])*V_i1*sind(teta2-teta1-atand.(imag(Y[2,1])/real(Y[2,1]))))+
#(abs(Y[2,3])*V_i3*sind.(teta2-teta3-atand(imag(Y[2,3])/real(Y[2,3])))))


# Jacobian for the partials of P2

dP2(V_i2,V_i3,teta2,teta3) = V_i2*((abs(Y[2,1])*V_i1*cos(teta2-teta1-atan.(imag(Y[2,1])/real(Y[2,1]))))+
(abs(Y[2,2])*V_i2*cos.(-atan(imag(Y[2,2])/real(Y[2,2]))))+(abs(Y[2,3])*V_i3*cos.(teta2-teta3-atan(imag(Y[2,3])/real(Y[2,3])))))

display("------------DP2--------")

dP2_dteta22(V_i2,V_i3,teta2,teta3) =(dP2.(V_i2,V_i3,teta2+h,teta3)-dP2.(V_i2,V_i3,teta2,teta3))/h
dP2_dV22(V_i2,V_i3,teta2,teta3)=(dP2.(V_i2+h,V_i3,teta2,teta3)-dP2.(V_i2,V_i3,teta2,teta3))/h
dP22_dteta3(V_i2,V_i3,teta2,teta3)=(dP2.(V_i2,V_i3,teta2,teta3+h)-dP2.(V_i2,V_i3,teta2,teta3))/h
dP22_dV3(V_i2,V_i3,teta2,teta3)=(dP2.(V_i2,V_i3+h,teta2,teta3)-dP2.(V_i2,V_i3,teta2,teta3))/h

display(dP2_dteta22(V_i2,V_i3,teta2,teta3))
display(dP2_dV22(V_i2,V_i3,teta2,teta3))
display(dP22_dteta3(V_i2,V_i3,teta2,teta3))
display(dP22_dV3(V_i2,V_i3,teta2,teta3))

#Jacobian for the partials of Q2

dQ2(V_i2,V_i3,teta2,teta3) = V_i2*((abs(Y[2,1])*V_i1*sin(teta2-teta1-atan.(imag(Y[2,1])/real(Y[2,1]))))+
(abs(Y[2,2])*V_i2*sin.(-atan(imag(Y[2,2])/real(Y[2,2]))))+(abs(Y[2,3])*V_i3*sin.(teta2-teta3-atan(imag(Y[2,3])/real(Y[2,3])))))

dQ2_dteta22(V_i2,V_i3,teta2,teta3)=(dQ2.(V_i2,V_i3,teta2+h,teta3)-dQ2.(V_i2,V_i3,teta2,teta3))/h
dQ2_dV22(V_i2,V_i3,teta2,teta3) =(dQ2.(V_i2+h,V_i3,teta2,teta3)-dQ2.(V_i2,V_i3,teta2,teta3))/h
dQ3_dteta3(V_i2,V_i3,teta2,teta3)=(dQ2.(V_i2,V_i3,teta2,teta3+h)-dQ2.(V_i2,V_i3,teta2,teta3))/h
dQ3_dV3(V_i2,V_i3,teta2,teta3)=(dQ2.(V_i2,V_i3+h,teta2,teta3)-dQ2.(V_i2,V_i3,teta2,teta3))/h

display("------------DQ2--------")
display(dQ2_dteta22(V_i2,V_i3,teta2,teta3))
display(dQ2_dV22(V_i2,V_i3,teta2,teta3))
display(dQ3_dteta3(V_i2,V_i3,teta2,teta3))
display(dQ3_dV3(V_i2,V_i3,teta2,teta3))


#Jacobian for the P3

k = 3
r = 1;
dP3(V_i2,V_i3,teta2,teta3) = V_i3*((abs(Y[k,r])*V_i1*cos(teta3-teta1-atan.(imag(Y[k,r])/real(Y[k,r]))))+
(abs(Y[k,r+1])*V_i2*cos.(teta3-teta2-atan(imag(Y[k,r+1])/real(Y[k,r+1]))))+(abs(Y[k,r+2])*V_i3*cos.(teta3-teta3-atan(imag(Y[k,r+2])/real(Y[k,r+2])))))

display("------------DP3--------")

dP33_dteta22(V_i2,V_i3,teta2,teta3) =(dP3.(V_i2,V_i3,teta2+h,teta3)-dP3.(V_i2,V_i3,teta2,teta3))/h
dP33_dV22(V_i2,V_i3,teta2,teta3) =(dP3.(V_i2+h,V_i3,teta2,teta3)-dP3.(V_i2,V_i3,teta2,teta3))/h
dP33_dteta3(V_i2,V_i3,teta2,teta3)=(dP3.(V_i2,V_i3,teta2,teta3+h)-dP3.(V_i2,V_i3,teta2,teta3))/h
dP33_dV3(V_i2,V_i3,teta2,teta3)=(dP3.(V_i2,V_i3+h,teta2,teta3)-dP3.(V_i2,V_i3,teta2,teta3))/h

display(dP33_dteta22(V_i2,V_i3,teta2,teta3))
display(dP33_dteta3(V_i2,V_i3,teta2,teta3))
display(dP33_dV22(V_i2,V_i3,teta2,teta3))
display(dP33_dV3(V_i2,V_i3,teta2,teta3))

# Jacobian for the Q

dQ3(V_i2,V_i3,teta2,teta3) = V_i3*((abs(Y[k,r])*V_i1*sin(teta3-teta1-atan.(imag(Y[k,r])/real(Y[k,r]))))+
(abs(Y[k,r+1])*V_i2*sin.(teta3-teta2-atan(imag(Y[k,r+1])/real(Y[k,r+1]))))+(abs(Y[k,r+2])*V_i3*sin.(teta3-teta3-atan(imag(Y[k,r+2])/real(Y[k,r+2])))))


dQ33_dteta22(V_i2,V_i3,teta2,teta3) =(dQ3.(V_i2,V_i3,teta2+h,teta3)-dQ3.(V_i2,V_i3,teta2,teta3))/h
dQ33_dV22(V_i2,V_i3,teta2,teta3) =(dQ3.(V_i2+h,V_i3,teta2,teta3)-dQ3.(V_i2,V_i3,teta2,teta3))/h
dQ33_dteta3(V_i2,V_i3,teta2,teta3)=(dQ3.(V_i2,V_i3,teta2,teta3+h)-dQ3.(V_i2,V_i3,teta2,teta3))/h
dQ33_dV3(V_i2,V_i3,teta2,teta3)=(dQ3.(V_i2,V_i3+h,teta2,teta3)-dQ3.(V_i2,V_i3,teta2,teta3))/h
display("------------DQ3--------")

display(dQ33_dteta22(V_i2,V_i3,teta2,teta3))
display(dQ33_dteta3(V_i2,V_i3,teta2,teta3))
display(dQ33_dV22(V_i2,V_i3,teta2,teta3))
display(dQ33_dV3(V_i2,V_i3,teta2,teta3))


J5(V_i2,V_i3,teta2,teta3) = [dP2_dteta22(V_i2,V_i3,teta2,teta3) dP22_dteta3(V_i2,V_i3,teta2,teta3) dP2_dV22(V_i2,V_i3,teta2,teta3) dP22_dV3(V_i2,V_i3,teta2,teta3);
dP33_dteta22(V_i2,V_i3,teta2,teta3) dP33_dteta3(V_i2,V_i3,teta2,teta3) dP33_dV22(V_i2,V_i3,teta2,teta3) dP33_dV3(V_i2,V_i3,teta2,teta3) ;
dQ2_dteta22(V_i2,V_i3,teta2,teta3) dP33_dV22(V_i2,V_i3,teta2,teta3) dQ2_dV22(V_i2,V_i3,teta2,teta3)  dQ3_dV3(V_i2,V_i3,teta2,teta3);
dQ33_dteta22(V_i2,V_i3,teta2,teta3) dQ33_dteta3(V_i2,V_i3,teta2,teta3) dQ33_dV22(V_i2,V_i3,teta2,teta3)  dQ33_dV3(V_i2,V_i3,teta2,teta3)
]


y_initial1(V_i2,V_i3,teta2,teta3) = [dP2.(V_i2,V_i3,teta2,teta3);dP3.(V_i2,V_i3,teta2,teta3);dQ2.(V_i2,V_i3,teta2,teta3);dQ3.(V_i2,V_i3,teta2,teta3)]

P2 =-2
P3 = 1
Q2 = -0.5
Q3 =
dely1 = [P2;P3;Q2]

y = [P2;P3;Q2]
delY = y_initial1(V_i2,V_i3,teta2,teta3)-y; # for my NR iterations


tol = 10^-2
counter = 0;
X_in = [teta2;teta3;V_i2;V_i3]
for i in 1:1

    X_t = ((J5.(X_in[3],X_in[4],X_in[1],X_in[2])\(delY)))
    X_next = X_in-X_t # book does it in an inverse way keep that in mind
    display(X_next[3])
    display("------------")
    display(X_next[2])
        if(abs(X_next[3]-X_in[3])<=tol)
            display(["Iterations : ",i])
            display(X_next[3]);
            break;
        end
    X_in =X_next;
    delY = y_initial1(V_i2,V_i3,teta2,teta3)-y; # for my NR iterations
end


J5_new(V_i2,V_i3,teta2,teta3) = [dP2_dteta22(V_i2,V_i3,teta2,teta3) dP22_dteta3(V_i2,V_i3,teta2,teta3) dP2_dV22(V_i2,V_i3,teta2,teta3) ;
dP33_dteta22(V_i2,V_i3,teta2,teta3) dP33_dteta3(V_i2,V_i3,teta2,teta3) dP33_dV22(V_i2,V_i3,teta2,teta3) ;
dQ2_dteta22(V_i2,V_i3,teta2,teta3) dP33_dV22(V_i2,V_i3,teta2,teta3) dQ2_dV22(V_i2,V_i3,teta2,teta3)  ;

]




Q3 = 0.241559124124






f1(x) = cos(x)
h = 10^-6
f2(x) = (f1(x+h)-f1.(x))*h.^-1





Pi = Vi_2*(Y)
