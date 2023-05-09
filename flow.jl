#
# value of Flow
F = 5:5:70;
#last value of deltaP
delP_last = 85;
#relationship
k = F[length(F)]^2/delP_last
#finding deltaP
deltaP = F.^2*(k^-1)

display(deltaP)

#4-20mA relationship

using SymPy

# for flow ->4-20mA
y =a-> 4+(20-4)/(F[length(F)]-F[1])*(a-F[1])

#for the pressure 4-20mA

y1 =a-> 4+(20-4)/(deltaP[length(deltaP)]-deltaP[1])*(a-deltaP[1])

fourtotwenty_flow = 0;
fourtotwenty_pressure = 0;
for i in 1:length(deltaP)

fourtotwenty_flow = [fourtotwenty_flow y.(F[i])]
fourtotwenty_pressure = [fourtotwenty_pressure y1.(deltaP[i])]


end
fourtotwenty_flow1 = fourtotwenty_flow[2]
fourtotwenty_pressure1 = fourtotwenty_pressure[2]
for i in 3:length(fourtotwenty_flow)
fourtotwenty_flow1 = [fourtotwenty_flow1 fourtotwenty_flow[i]]
fourtotwenty_pressure1 = [fourtotwenty_pressure1 fourtotwenty_pressure[i]]

end


fourtotwenty_flow1 = vec(fourtotwenty_flow1);
fourtotwenty_pressure1 = vec(fourtotwenty_pressure1);

display(fourtotwenty_flow)
display(fourtotwenty_pressure)
