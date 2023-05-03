%random example for comparing with linprog
A = [50 30 30 1 0 0; 2 3 2 0 1 0; 1 1 1 0 0 1];
b = [2000; 70; 30];
c =[20000 15000 16000 0 0 0]';
Beta = [4 5 6];
v = [1 2 3];

%running the file simplex which is our implementation of the simplex method
simplex(A, b, c, Beta, v);

%linprog implementation
A = [50 30 30; 2 3 2; 1 1 1];
c = [-20000 -15000 -16000];
b = [2000 70 30];
Aeq = [];
Beq = [];
lb = [0 0 0];
ub = [inf inf inf];
[X, Z] = linprog(c, A, b, Aeq, Beq, lb, ub);
cost_linprog = Z