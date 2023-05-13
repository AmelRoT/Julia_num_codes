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



#f(x,y) = 16-x^2-2*y^2
f = f1; 
V = 0;  #volume 

for j in 1:length(y_f)

    for i = 1:length(x_f) 

        V = V+f(x_f[i],y_f[j])*delta_X*delta_Y;
    end 
end 

    display(V); 

end 



f1(x,y) = (1-x^2)^(1/2)

#f1(x,y) = 16-x^2-2*y^2
doubleIntegralWithEdgePoints(f1,500,500,-1,1,-2,2)

















