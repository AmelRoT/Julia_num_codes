# taylor series for exp (x) centered at 0

x = 1;
a = 0;
m  = 1;
tay = 1;
n = 1;
t = 1
sum = 1
c = 30
for i in 1:c

t = x*t
n = n*(m);
m = m+1;
sum = sum+t/n
a = [a sum]
if(m==c)
    break;
end

end
