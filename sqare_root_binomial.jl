


#---------------------binioanl sqrt-------------------------#

# (1+x)^n = sum (1+nx+n(n-1)/2!+...+)

result = 1; 
n = 1/2; 

result= x-> 1+(n)*x+(n*(n-1)/2)*x^2+(n*(n-1)*(n-2)/6)*x^3+(n*(n-1)*(n-2)*(n-3)/24)*x^4#(n*(n-1)*(n-2)*(n-3)*(n-4)/120)*x^5+(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)/(120*5))*x^6
n = 1/2
r(x) = 1+(n)*x+(n*(n-1)/2)*x^2+(n*(n-1)*(n-2)/6)*x^3+(n*(n-1)*(n-2)*(n-3)/24)*x^4+(n*(n-1)*(n-2)*(n-3)*(n-4)/120)*x^5+(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)/(120*5))*x^6



n1 = 1; 
x = 3
n = 1/2
result = 1
for k  = 1:1:110

    for i = 0:1:(k-1)
        n1 = n1*(n-i)
    end 
    display(n1)
    display(i)
    result = result+(n1/(factorial(big(k))))*x^(k)
    n1 = 1; 


end 

x = 1
n1 = 1000
for n = 1:30

    x = (n1+x^2)/(2*x)

end 
display(x)


x = 1
n1 = 27
for n = 1:30

    x = (n1+x^3)/(2*x^2)

end 
display(x)


exp(log(9)/2)

5^-3


x = 2
e1 = 1; 
m = 1;
x1 = 1 
for i = 1:1:30
m = m*i; #using it like a factorial 

        x1= x1*x; 
 
    e1 = e1+x1/m
    display(e1)
end 

exp(2)


x1 = 1; 
x = 0.3; 
l1 = 0; 

for i = 1:1:30

    x1 = x1*x; 
    
    l1 = l1-x1/i; 
    display(l1)

end 


x1 = 1; 
x = 2; 
l1 = 0; 

for i = 1:1:20

    x1 = x1*x; 
    if(i%2 ==0)
        x1 = -x1; 
    end 
    if(i%2 != 0)
        x1 = abs(x1); 
    end 
    l1 = l1+x1/i; 
    display(l1)

end 

l1 = 0; 
x = 0.23
for n = 1:1:50

    l1 = l1+(((x-1)/(x+1)).^(2n-1))/(2n-1)


end 
l1 = 2*l1

7.389056098930649-exp(2)


log(5)

x = 5
f1 = n->((x-1)/(x+1)).^(2n-1)
log(5)-1.6094379124330567

1.6094379124341-log(5)

2^-2.5