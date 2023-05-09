using LinearAlgebra;

A = [4 12 8 4;1 7 18 9;2 9 20 20;3 11 15 14];
B = [4 ;-6;7;10];

L = 1;
for i in 1:length(A[1,:])

for j in 1:length(A[1,:])

    if (i==j)
    L = [L 1];
    end


    if(i!=j)
    L = [L 0];

    end

end

end

L1 = L[2];
for i in 3:length(L)

L1 = [L1 L[i]]

end
L2 = zeros(length(A[1,:]),length(A[1,:]))

for i in 1:length(L1)

L2[i] =L1[i]

end

L = L2;

U = A;
m1 = 0
for iter = 1:length(A[1,:])-1

m1 = m1+iter
end
mult = 1;
m2 = 2
m4 = m2;
m3 = 1
m5 = m3;
var = 0
n = 1;

for counter1 in 1:(length(A[1,:])-1)

for i in m2:(length(A[1,:]))


for j in m3:(length(A[1,:]))
if(m3==j)
mult = U[i,j]/(U[i-n,j])
L[i,j] = mult;
var = [var mult]
end

U[i,j] = U[i,j] - mult*U[i-n,j]

end
m2 = m2+1;
n = n+1;
end
m2 = m4+1;
m4 = m4+1;
m3 = m5+1;
m5 = m5+1;
n = 1;
end

var1 = var[2]

for i in 3:length(var)
var1 = [var1 var[i]]

end

# We want to find Ax = B -> A = L*U
# L*U*x= B -> U*x=Y -> L*Y = B

B = float(B);
y = zeros(length(B));
y[1] = B[1]/L[1,1];
v = 0;
for i in 1:length(B);

for j in 1:length(B)
    if(i!=j)
    v = L[i,j]*y[j]+v
    end
end
if(i!=1)
s = (B[i]-v)/(L[i,i])
v = 0;
y[i]=s;
end
end

# we got Y -> U*X =Y

x = zeros(length(y));
x[length(y)] = y[length(y)]/U[length(y),length(y)];
v = 0;
for i in 1:length(y);

for j in 1:length(y)
    if(length(y)!=i)
    v = U[length(y)-i,j]*x[j]+v
    end
end
if(i!=length(y))
s = (y[length(y)-i]-v)/(U[length(y)-i,length(y)-i])
v = 0;
x[length(y)-i]=s;
end
end

display(x)
