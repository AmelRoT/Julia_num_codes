#------------------------Testing Square Root---------#

n = 300
result = 0;
i = 0
while(true)

    if((i*i)>n)
        display(i)
        result = i-1;
        break;
    end
    i = i+1;
end

while(true)

    result = result+0.000001;


    if(result*result>n)
        result = result-0.000001;
        break;
    end

end

display(result)


    #------------------Based on Bionmial Expansion ---------#

function sqroot(x)

    n = 10;
end
    fact(10)

    
