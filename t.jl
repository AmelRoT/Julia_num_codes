using Plots
function splin(x,y)
    # z=[z_1=0, z_2, ..., z_(n-1), z_n=0]
    n = length(x)
    A = zeros(n-2,n-2)
    h = [x[i+1]-x[i] for i in 1:1:n-1]
    b = [(y[i+2]-y[i+1])/h[i+1]-(y[i+1]-y[i])/h[i] for i in 1:1:n-2]
    b = 6*b
    for i in 2:1:n-3
        A[i,i-1] = h[i]
        A[i,i]   = 2*(h[i]+h[i+1])
        A[i,i+1] = h[i+1]
    end
    A[1,1] = 2*(h[1]+h[2]); A[1,2] = h[2]
    A[n-2,n-3] = h[n-2]; A[n-2,n-2] = 2*(h[n-2]+h[n-1])
    z =  A\b
    z = [0;z;0]
    display(A)
    display(b)
    q(j) = xx ->((z[j]/(6*h[j]))*(x[j+1]-xx)^3 + (z[j+1]/(6*h[j]))*(xx-x[j])^3 +
                (y[j]/h[j]-(z[j]*h[j])/6)*(x[j+1]-xx) +
                (y[j+1]/h[j]-(z[j+1]*h[j])/6)*(xx-x[j]))

    scatter(x,y,xlim=[x[1],x[end]],label="(x,y)")
    for k in 1:1:n-1
        xq = x[k]:0.01:x[k+1]
        yq = q(k).(xq)
        Plots.display(plot!(xq,yq,label="q($k)",framestyle=:origin))
        sleep(0.1)
    end
end
#############################################################################


x =[8 11 15 18 22]
y = [5 9 10 8 7]
x= [4.6  4.8  5.0  5.2 5.4 5.6 5.8 6.0]
y = [32.11 9.00 -3.52 -9.23 -11.58 -12.01 -11.24 -10.12]
x = 0:1:6
y = sin.(x)
splin(x,y)
xx = 0:0.01:6
yy = sin.(xx)
plot!(xx,yy,label="sin(x)",lw =2,ls=:dash)
