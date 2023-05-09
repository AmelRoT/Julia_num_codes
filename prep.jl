
#-------------------Flow from DiffP-----------------------
#flow parameters
yi = 10;
yf = 30;
#differential pressure parameters
xf = 10;
xi = 1;

inputV = xi:1:xf; #intial values for the DiffP

y =x-> yi.^2+ (yf.^2-yi.^2)/(xf-xi)*(x-xi)
ynew = y.(inputV)

q = sqrt.(ynew) #flow is sqrt of y

using Plots
plotly()

scatter(inputV,q)







#------------From DiffP to Flow ------------------------------------
#diff P  parameters
yi = 50;
yf = 100;
# flow parameters
xf = 30;
xi = 0;

inputV =0:3:30  #for flow

y =x-> yi+ (yf-yi)/(xf.^2-xi.^2)*(x.^2-xi.^2)
ynew = y.(inputV) #differential pressure

using Plots
plotly()

scatter(inputV,ynew)

#-------------------4-20mA diffP------------------

xnew = ynew;
#Current parameters
yi1 = 4;
yf1 = 20;

#Diff parameters
xf1 = 100;
xi1 = 50;

y1 =x-> yi1+(yf1-yi1)/(xf1-xi1)*(x-xi1)
y_current = y1.(xnew)
display(y_current)

scatter(xnew,y_current)





#---------------------4-20mA flow---------------------------
xnew = inputV;
#Current parameters
yi1 = 4;
yf1 = 20;

#flow parameters
xf1 = 30;
xi1 = 0;

y1 =x-> yi1+(yf1-yi1)/(xf1-xi1)*(x-xi1)
y_current = y1.(xnew)
display(y_current)

scatter(xnew,y_current)
