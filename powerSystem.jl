#-----------------Gauss Siedel mehtod and other -------------#

A1 = [10 5; 2 9]
B1 = [6;3]

y = inv(A1)*B1;
display(y[1]);
display(y[2])

#---------------------Next method ------------------#
A1 = [2 3 -1;-4 6 8;10 12 14]
B1 = [5;7;9]
y = inv(A1)*B1;
display(y[1])
display(y[2])
display(y[3])

##
using SymPy;


Zm= symbols("Zm");
Zn = symbols("Zn");
Zs = symbols("Zs");
Z1 = Zs+Zn;
Z2 = Zm+Zn
Z = [Z1 Z2 Z2;Z2 Z1 Z2; Z2 Z2 Z1]

A = [1 1 1; 1 -0.5000 - 0.8660im -0.5000 + 0.8660im;1 -0.5000 + 0.8660im -0.5000 - 0.8660im ]

ZZ = inv(A)*Z*A
Zz = simplify(ZZ)
