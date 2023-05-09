
#----------------------------Monte Carlo approximation ---------------------------

#------------------------------Generating random numbers ------------------------
using Plots
plotly()

x = rand(-1:0.0001:1)

y =  rand(-1:0.0001:1);
res = 0;
n = 8000;
xu = 0;
yu = 0;


for i in 1:n
    x = rand(-1:0.0001:1)
    y =  rand(-1:0.0001:1);
if(x^2+y^2<=1)
xu = [xu x]
res = res+1;
yu = [yu y]
end

end

A_circle = 4*res/n

x1 = -1:0.001:1

y1 =1
y2 = -1

y3 = -1:0.001:1;
x3 = -1;
x4 = 1
for i in 1:length(x1)-1

y1 = [y1 1]
y2 = [y2 -1]
x3 = [x3 -1]
x4 = [x4 1]
end

#circle
y5 = (1-x1[1].^2).^(1/2)
for i in 2:length(x1)

y5 =[y5 (1-x1[i].^2).^(1/2)]


end


31.4*6

r = 1;
theta = 0:pi/40:2*pi
xnew= r*cos.(theta)
ynew = r*sin.(theta)

plot(x1,vec(y1),linewidth="3",color="blue")
plot!(x1,vec(y2),linewidth="3",color="blue")
plot!(vec(x3),y3,linewidth="3",color="blue")
plot!(vec(x4),y3,linewidth="3",color="blue")
plot!(xnew,ynew,linewidth="3",color="red")

scatter!(xu,yu)
scatter!(xu,yu)

#what is the primary reason for the breakpoint ----example driven part ->

for i in 1:n
    x = rand(-1:0.0001:1)
    y =  rand(-1:0.0001:1);
if(x^2+y^2<=1)
xu = [xu x]
res = res+1;
yu = [yu y]
end


for i in 1:n

x = rand(-1:0.01:2)

if(x^2+x*y.^2 = 100)
scatter!(x,y)
plot!()

end
end


#--------------------------------------------------

x = rand(0:0.0001:3)
y = rand(0:0.0001:2)

res = 0;
n = 50000;
xu = 0;
yu = 0;


for i in 1:n

    x = rand(0:0.0001:3)
    y = rand(0:0.0001:9)
if(x.^2>=y)
xu = [xu x]
res = res+1;
yu = [yu y]
end

end

A_circle = 27*res/n

x2= -2:0.0001:2
y2 = x2.^2
plot(x2,y2,linewidth="3",color="blue")



x1 = 0:0.001:4

y1 = 0.0
y4 = 4.0
y3 = 0:0.001:4;
x3  = 0.0;
x4 = 2.0;
for i in 2:length(x1)

y1 = [y1 0.0]
y4 = [y4 4.0]

end

for i in 2:length(y3)


    x3 = [x3 0.0]
    x4 = [x4 2.0]

end
plot!(x1,vec(y1),linewidth="3",color="red")
plot!(x1,vec(y4),linewidth="3",color="red")
plot!(vec(x3),y3,linewidth="3",color="red")
plot!(vec(x4),y3,linewidth="3",color="red")

scatter!(xu,yu)


#-------------------------------------------------

x = rand(-3:0.0001:3)
y = rand(0:0.0001:1)

res = 0;
n = 15000;
xu = 0;
yu = 0;


for i in 1:n

    x = rand(-3:0.00001:3)
    y = rand(0:0.00001:1)

if(exp.(-x.^2)>=y)
xu = [xu x]
res = res+1;
yu = [yu y]
end

end

A_circle = 6*res/n

x6 = -3:0.001:3
y6 = exp.(-x6.^2)

plot(x6,y6,linewidth="3",color="blue")
x1 = -3:0.001:3

y1 = 0.0
y4 = 1.0
y3 = 0:0.001:1;
x3  = 3.0;
x4 = -3.0;
for i in 2:length(x1)

y1 = [y1 0.0]
y4 = [y4 1.0]

end

for i in 2:length(y3)


    x3 = [x3 3.0]
    x4 = [x4 -3.0]

end
plot!(x1,vec(y1),linewidth="3",color="red")
plot!(x1,vec(y4),linewidth="3",color="red")
plot!(vec(x3),y3,linewidth="3",color="red")
plot!(vec(x4),y3,linewidth="3",color="red")
scatter!(xu,yu)

#-------------------------------------------------------
using Plots
plotly()
x = [1 2 3 5 6]
y = [1 2 3 5 6]
z = [1 2 3 5 6]


plot(x,y,z)



x = rand(-1:0.001:1)
y = rand(-1:0.001:1)
z = rand(-1:0.001:1)


res = 0;
n = 30;
xu = 0;
yu = 0;
zu = 0;

for i in 1:n

    x = rand(-1:0.001:1)
    y = rand(-1:0.001:1)
    z = rand(-1:0.001:1)

if(x.^4+y.^4+z.^4<=1)
xu = [xu x]
res = res+1;
yu = [yu y]
zu = [zu z]
end

end

V_n = 8*res/n

x1 = -1:0.0001:1
y1 = 1.0
y2 =-1.0
y3 = -1:0.0001:1
z1 = -1;
z2 = 1;
for i in 2:length(x1)

y1 = [y1 1.0]
y2 = [y2 -1.0]
z1 = [z1 -1]
z2 = [z2 1]
end

plot(x1,vec(y1),vec(z1),linewidth="3",color="blue")
plot!(x1,vec(y2),vec(z1),linewidth="3",color="blue")
plot!(vec(y1),x1,vec(z1),linewidth="3",color="blue")
plot!(vec(y2),x1,vec(z1),linewidth="3",color="blue")



plot!(x1,vec(y1),vec(z2),linewidth="3",color="blue")
plot!(x1,vec(y2),vec(z2),linewidth="3",color="blue")
plot!(vec(y1),x1,vec(z2),linewidth="3",color="blue")
plot!(vec(y2),x1,vec(z2),linewidth="3",color="blue")




scatter(xu,yu,zu)

scatter(x1,y3)



#------------------------------------------------------------

res = 0;
n = 8000;
xu = 0;
yu = 0;
zu = 0;

for i in 1:n

    x = rand(-1:0.0001:1)
    y = rand(-1:0.0001:0)
    z = rand(0:0.0001:1)

if(x.^2+y^2<=1 && z<=-y )
xu = [xu x]
res = res+1;
yu = [yu y]
zu = [zu z]
end

end

V_n = 2*res/n

x1 = -1:0.0001:1
y1 = 1.0
y2 =-1.0
y3 = -1:0.0001:1
z1 = -1;
z2 = 1;
for i in 2:length(x1)

y1 = [y1 1.0]
y2 = [y2 -1.0]
z1 = [z1 -1]
z2 = [z2 1]
end

plot(x1,vec(y1),vec(z1),linewidth="3",color="blue")
plot!(x1,vec(y2),vec(z1),linewidth="3",color="blue")
plot!(vec(y1),x1,vec(z1),linewidth="3",color="blue")
plot!(vec(y2),x1,vec(z1),linewidth="3",color="blue")



plot!(x1,vec(y1),vec(z2),linewidth="3",color="blue")
plot!(x1,vec(y2),vec(z2),linewidth="3",color="blue")
plot!(vec(y1),x1,vec(z2),linewidth="3",color="blue")
plot!(vec(y2),x1,vec(z2),linewidth="3",color="blue")




scatter!(xu,yu,zu)

scatter(x1,y3)


#----------------------------------------------------------------------
using Plots
V1 = 0
res = 0;
n = 100000;
xu = 0;
yu = 0;
zu = 0;

for i in 1:n

    x = rand(-1:0.0001:1)
    y = rand(-1:0.0001:1)
    z = rand(-1:0.0001:1)

if(abs(x).^0.5+abs(y)^0.5+abs(z).^0.5<=1)
xu = [xu x]
res = res+1;
yu = [yu y]
zu = [zu z]
end

end


V_n = 8*res/n

V1 = [V1 V_n]


scatter(xu,yu,zu)
