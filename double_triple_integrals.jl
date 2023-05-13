# double integrals with square limits 

 
# constructing the square 

#     __________________
#  d |                  |
# ym |                  |
#... |                  |
# y1 |                  |
#  c | _ _ _ _ _ _ _ _ _| 
#    a    x1 ...  x_n   b


function doubleIntegralWithEdgePoints(f1,n,m,a,b,c,d)


    delta_X = (b-a)/n
    delta_Y = (d-c)/m
    x_f = a+delta_X:delta_X:b; 
    y_f = c+delta_Y:delta_Y:d; 

    f = f1; # function 
    V = 0;  # volume 

    for j in 1:length(y_f)

        for i = 1:length(x_f) 

            V = V+f(x_f[i],y_f[j])*delta_X*delta_Y;
        end 
    end 

        display(V); 

end 



f1(x,y) = (1-x^2)^(1/2)
f1(x,y) = x*sin(y)


f1(x,y) = (2*x+y)^8
f1(x,y) = (x*exp(x))/(y)
f1(x,y) = x*(sin(y)).^2

#f1(x,y) = 16-x^2-2*y^2
doubleIntegralWithEdgePoints(f1,500,500,-1,1,-2,2)
doubleIntegralWithEdgePoints(f1,500,500,0,2,0,pi/2)
doubleIntegralWithEdgePoints(f1,2*10^3,2*10^3,0,1,0,2)
doubleIntegralWithEdgePoints(f1,2*10^3,2*10^3,0,1,1,2)
doubleIntegralWithEdgePoints(f1,2*10^3,2*10^3,0,2,0,pi)




function doubleIntegralTypeI(f1,n,m,a,b,c,d,y_up, y_do)


    delta_X = (b-a)/n
    delta_Y = (d-c)/m
    x_f = a+delta_X:delta_X:b; 
    y_f = c+delta_Y:delta_Y:d; 

    y_upper = y_up; 
    y_down = y_do; 

    f = f1; # function 
    V = 0;  # volume 

    for j in 1:length(y_f)

        for i = 1:length(x_f) 
            if((y_f[j]<=y_upper(x_f[i])) && (y_f[j]>= y_down(x_f[i])))
               V = V+f(x_f[i],y_f[j])*delta_X*delta_Y;
            end 
        end 
    end 

        display(V); 

end 

y1(x) = 1+x^2
y2(x) = 2*x^2
f1(x,y) = x+2*y

doubleIntegralTypeI(f1,2*10^3,2*10^3,-1,1,0,10,y1,y2)

















