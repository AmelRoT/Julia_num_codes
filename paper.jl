using Plots
plotly()


# Eulter method

a = 0
b = 0.6
n = 100

h = (b-a)/n; # step size

# setting up the x ->
f(x,y) = 3*x*exp.(-y)

function EulerMethod(f,a,b,n,yinit,xi) #f (function,limits,number of steps,yinit)
h = (b-a)/n
#xi = a:h:b

yi = yinit; #intial cond
s = yi

for i = 2:length(xi)

 yf = yi+h*f.(xi[i],yi);
 yi = yf
 s = [s yi]

end
return s
end
s1 = EulerMethod(f,a,b,n,2)
scatter!(xi,vec(s1))

#-------------------Backward Euler method ----------------------------#

#-------------------------------------Newton method modified -------------------------#
function NewtonMethod(f,xo,j1)
        h1 = 10^-5
        j = j1;
        f1(x,x2) = (f.(x+h1,x2)-f(x,x2))/h1; # approximation of derivative
        for i in 1:300
                x_n1 = xo-f.(xo,xi[j])*(f1.(xo,xi[j])).^-1
                xo = x_n1

                    display(x_n1)
                    if(i==100)
                            return x_n1;
                    end
        end
end

#-------------------------------------------------------------------------#




a = 0;
b = 0.6
n = 100;

f2(x,y) = 3*x*exp.(-y)


function BackwardEuler(f2,a,b,n,yinit)
        h = (b-a)/n
        xi = a:h:b
        yi =yinit
        y(yf,x) = -yf+yi+h*f2.(x,yf)

        s = yi;
        for i = 2:length(xi)
        Newt = NewtonMethod(y,yi,i)
        s = [s Newt];
        yi = Newt;

end
       return s;
end


scatter(xi,vec(s))

#------------------------------#

a = 0
b = 2

n = 100;
h = (b-a)/n
xi = a:h:b

m = 10; #const
alfa = 2 # const
g = -9.81
k = 30;

U_c = -1/m*ones((n+1),(n+1))
Q = alfa*ones((n+1),(n+1))

f5(v,lambda) = v+k/m*lambda

f3(lambda,v) = g+(1/(m*m*alfa))*lambda-(k/m)*v
v_init = 2
res1 = BackwardEuler(f3,a,b,n,v_init);
scatter(xi,vec(res1))

f4(res1,lambda) = res1+k/m*lambda
lam_init = 0
res2 = EulerMethod(f4,a,b,n,lam_init,res1)
scatter(xi,vec(res2))

#-----------------------------Paper system --------------------------#

a= 0;
b = 2;
N = 100;

#initial cond
Vin = 2 #Initial Vin
g = -9.81
deltaT = (b-a)/N;
m = 10;
alfa = 1
k = 10


c11 = -Vin;
c13 = 0;
        #-------------------Generating C11 and C13-----------#

for i = 2:N

        c11 =[c11 0]
        c13 = [c13 0]
end
        c13 = [c13 0]
        c11 = vec(c11)
        c13 = vec(c13)

        # ------------------Generating C12-----------#

c12 = g+(Vin/deltaT)
for i = 2:N

        c12 =[c12 g]

end
        c12 = vec(c12)


#-----------------------------Generating R matrix ------------------#

R = zeros((N),(N))

for i = 1:(N)

for j = 1:(N)

        if(i==j && i!=N)
                R[(i+1),j] = 1
        end

end
end

#-----------------------------Generating U1 and U2 matrix ------------------#

U1 = zeros((N+1),(N))

for i = 1:(N)

for j = 1:(N)

        if(j==(i))
                U1[i,j] = -1/m
        end

end
end

U2 = zeros((N),(N+1))

for i = 1:(N)

for j = 1:(N)

        if(j==(i))
                U2[i,(j+1)] = -1/m
        end

end
end

using LinearAlgebra
U1 = zeros(N+1,N)
 U1[1:N,:] = (-1/m)*Matrix{Float64}(I,N,N)

 U2 = zeros(N,N+1)

U2[1:N,(2:N+1)] = (-1/m)*Matrix{Float64}(I,N,N)







#-----------------------------Generating Q matrix ------------------#

Q = zeros((N+1),(N+1))

for i = 1:(N+1)

for j = 1:(N+1)

        if(i==j)
                Q[i,j] = alfa
        end

end
end

#-----------------------------Generating A and A^T matrix ------------------#

A = zeros((N),(N))


for i = 1:(N)


        for j = 1:(N)

                if(j==i)
                        A[i,j] = ((k/m)+1/(deltaT))
                end

                if(j==(i+1))

                        A[i,j] = (-1/(deltaT))

                end


        end

end

 A_t = transpose(A)

A_total = [R A zeros(N,N+1);A' zeros(N,N) U2;zeros((N+1),N) U1 Q]

c_total = [c11;c12;c13]

X = inv(A_total)*c_total

v = X[1];

for i = 2:(N)


        v = [v X[i]]


end
        v = vec(v)

lambda = X[N+1];

for i = 2:(N)


        lambda = [lambda X[N+i]]


end
        lambda = vec(lambda)

u = X[N*2+1];

for i = 2:(N+1)


        u = [u X[(2*N)+i]]


end

        u = vec(u)


Anew = U1*inv(A)*R*inv(A)*U2+Q
B = c13+U1*inv(A)*R*inv(A')*c12-U1*inv(A)*c11

u1 =Anew\B


#------------Conjugate Gradient Method ------------------#

A = [1 23 5;3 4 31; 123 -2 100;]
b = [110;20.5; 30.123]


A1 = A'*A
b1 = (A')*b
x_n = [1;2;100]


A = A_total
b = c_total
x_n = c_total

A1 = A_total'*A_total
b1 = A_total'*c_total
x_n = b1*2;

function conjugateGradientMethod(A,b,x_n)
r = b-A*x_n
p = r

        for i = 1:length(b)
                alfa = (r'*r)/(p'*A*p)
                x_new = x_n+alfa*p
                r_new = r-alfa*A*p
                beta = (r_new'*r_new)/(r'*r)

                #display(x_new)
                if(sqrt((r'*r))<10^-10)

                        display(["I : " i])
                        break;
                end
                p_new = r_new+beta*p
                p = p_new
                r = r_new
                x_n = x_new
        end
        return x_n
end


#-------------------------GMRES-----------------------#
