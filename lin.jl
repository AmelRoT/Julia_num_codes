A =[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0]

    B= [0 0 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1]



A = [2 -3 1 2 4;3 -1 0 3 2;7 0 -1 7 2;4 1 -1 4 0;1 2 -1 1 -2]

A = [1 1 0 1 3 2;1 3 -1  0 1 1;2 2 0 2 6 4;1 -1 2 -3 0 3]

A= [-1 0 3 2 8;0 -1 7 2 2;1 -1 4 0 4;2 -1 1 -2 -6]

A1 = A'

A = [0 -1 2 2;1 -1 0 4;1 -2 2 6]

A = [3 2 2;1 -2 -10;-1 -2 -6]

A = [2 -2 8;1 -1 4;3 -3 12]

A = [2 -2 8; 1 -1 4;3 -3 2]

A = [2 2 0;1 -1 0;3 -3 1]

u = [0 1 0 -2 -2]
v =[0 2 -2 1 0]


A = [0 5 1;-1 2 1;2 1 -1;-1 3 3]

A = [5 2 -2 ; 1 -2 0;0 2 -3; 1 2 3]


A = [1 0 2;-2 5 -4;4 4 4]

A = [1 1 0;0 -2 5;2 0 5]


A = [1 1 -4;2 2 6;1 4 4]

A = [1 3 0 2;0 1 1 0; 0 0 1 -1;1 -1 0 1]


A = [2 1 4 2;0 3 2 8;2 4 6 10;2 7 8 18]

A = [-1 4 2 5 -6;-3 1 4 -2 7;3 -2 1 4 -5;1 -1 0 2 4]


A = [4 -5 2 -1;-3 1 4 -2;1 4 6 -3]

u = [4 -5 2 -1]
v = [-3 1 4 -2]
g = [1 4 6 -3]

A = [2 -2 8;1 -1 4;3 -3 12]

A = [2 -2 8;1 -1 4;3 -3 2]

A = [0 0 0;0 0 0;0 0 0]

A = [2 2 0;1 -1 0; 3 -3 1]


A = [1 -1 2;-3 3 -6]


A = [1 3 0 1;0 1 2 -1;1 4 2 0]

A = [1 2 -1 0 ;2 1 0 0;3 3 5 0]

A = [1 -1 2;-3 3 -6]


A = [1 3 0 1;0 1 2 -1;1 4 2 0]

A = [1 2 1;1 1 0; 0 0 0]

A = [1 3 0;-1 1 1; 0 4 1]

A = [1 2 -1 0;0 -1 1 4;1 1 0 4]

A = [-1 4 5 2 -1;-3 1 0 -2 1; 1 2 5 6 1;1 -1 0 2 1]

A = [1 0 2;-1 1 0]

A = [1 -1 2 4 ;2 0 1 1]

A = [1 4 5 6 9;3 -2 1 4 -1;-1 0 -1 -2 -1;2 3 5 7 8]

A = [1 -1 2 0; 1 0 3 1; 3 -2 7 1;0 1 0 1;]

A = [1 -1 2 0;1 0 3 1;3 -2 7 1;0 1 0 1]

using LinearAlgebra
using RowEchelon


A = [-6 6;-3 3]
x = inv(A)*[0;0]

A=  [4-l 1 0 ;3 4-l 1 ;0 1 4-l]
l = symbols("l")

d = det(A)

B = l->[-3-l -3 -1;0 1-l 3;0 0 2-l]

B = l-> [4-l 1 0 ;3 4-l 1 ;0 1 4-l]

A = [3-l 2 ;1 4-l]
A = [5-l 0;1 5-l]


A=  [1-l 2 -1;0 3-l -2 ;0 0 -1-l]
B = l->[1-l 2 -1;0 3-l -2 ;0 0 -1-l]

M = [-1 1 3;2 0 5;7 -3 1]

d = det(M);

M = [2 1 3 -9;1 -1 2 -7;4 3 5 -15]

M= [0 -1 2 -1;5 2 1 3;1 1 -1 3]

M = [5 -1 0 1;2 2 2 2;-2 0 -3 3]

M = [5 2 -2;-1 2 0;0 2 -3]

M = [5 2 -2;-1 2 0; 1 2 3]

M = [0 5 1;-1 2 1; -1 3 3]

M = [2 3 1;-1 1 0; 1 4 1]

det(M)

A = [1 2 1;-1 0 2;0 2 3]
a = rref(A);

A = [1 -1 2; -3 3 -6]

A = [1 3 0 1;0 1 2 -1; 1 4 2 0]

A = [1 2 -1;2 1 0;3 3 5]

A = [1 3 0;-1 1 1;0 4 1]

A = [3 1 1 1;5 -1 1 -1]

A = [2 1 4 2;0 3 2 8;2 4 6 10; 2 7 8 18;]

A = [2 -3 1 2 4;3 -1 0 3 2;7 0 -1 7 2; 4 1 -1 4 0;1 2 -1 1 -2]

B = [10 -1 33;5 -4 14;7 -61 2 ]

A = [2 -2 4 4;2 -3 1 1;-2 4 2 2;0 1 3 3]

A = [1 0 1 3 2;3 -1 0 1 1;2 1 -3 1 4;-1 2 -3 0 3]

A = [1 0 0 ;0 cosd.(30) sind.(30); 0 sind.(30) cosd.(30)]
using SymPy;
l = symbols("l")
A = [3-l 2;1 4-l]

sol = solve(det(A))

B = l->[0-l 0 3;1 0-l -1;0 1 3-l]


A = [4-l 1 0;3 4-l 1;0 1 4-l]
B = l->[4-l 1 0;3 4-l 1;0 1 4-l]
F = B.(2)

A = [0-l 0 3;1 0-l -1;0 1 3-l]

A = [1-l 2 -1;0 3-l -2;0 0 -1-l]

 B= l ->  [1-l 0 -2;0 0-l 0;-2 0 4-l]

F = B(-1)

A = [1-l 4 ; 2 3-l]
F = B(5)

A = [1-l 0 -2;0 0-l 0;-2 0 4-l]

F = rref(B.(2))

A  = [2-l 1 0; 0 2-l 3; 0 0 5-l]

B =l-> [2-l 1 0; 0 2-l 3; 0 0 5-l]

A = [-1-l 9;-2 8-l]

sol = solve(det(A))
display(sol)

B =l-> [-1-l 9;-2 8-l]

A = [4-l 1 0;3 4-l 1;0 1 4-l]

det([2 1 -1;1 0 1;3 1 2])

A = [1/2 -sqrt(3)/2;sqrt(3)/2 1/2]
B = [sqrt(2)/2 -sqrt(2)/2; sqrt(2)/2 sqrt(2)/2]

F = [-1 0 ; 0 1]

A = [1 2 0 1;-1 0 1 1;1 4 1 3]

A = [1 2 0 1;1 0 1 1;0 2 -1 1;2 4 0 3]

A= [1 0 1;1 1 -1;3 1 2]

A = [1 -1 2 0;1 0 3 1; 3 -2 7 1;0 1 1 1 ]

A = [1/2 -sqrt(3)/2; sqrt(3)/2 1/2]

B = [1 0 ; 0 -1]

A = [-1-l 6;-3 8-l]

B = l->[-1-l 6;-3 8-l]



A = [-4 7;7 -2]
