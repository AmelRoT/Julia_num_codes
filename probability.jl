# using Julia for probability
x1 = 1;
x2 = 2;
x31 = 3;
x4 = 4;
for i in 1:98
        x3 = rand([1,2,3,4]);
        if(x3==1)
                x1 = [x1 x3]
        end
         if(x3 ==2)
         x2 = [x2 x3]
        end
        if(x3 ==3)
        x31 = [x31 x3]
       end
       if(x3 ==4)
       x4 = [x4 x3]
      end
end
