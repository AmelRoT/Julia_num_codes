# double integrals with square limits 

 
# constructing the square 

#     __________________
#  d |                  |
# ym |                  |
#... |                  |
# y1 |                  |
#  c | _ _ _ _ _ _ _ _ _| 
#    a    x1 ...  x_n   b


n = 500
m = 500

a = 0 
b = 2

c = 0 
d = 2

delta_X = (b-a)/n
delta_Y = (d-c)/m
x_f = delta_X:delta_X:b; 
y_f = delta_Y:delta_Y:d; 



f(x,y) = 16-x^2-2*y^2
V = 0;  #volume 

for j in 1:length(y_f)

    for i = 1:length(x_f) 

        V = V+f(x_f[i],y_f[j])*delta_X*delta_Y;
    end 
end 
























