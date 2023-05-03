%random example for comparing with linprog
A = [1 2 2 4 1 0 0 0; 50 75 100 10 0 1 0 0;3 1 0 0 0 0 1 0; 8 3 1 1 0 0 0 1];
b = [250; 8000; 4500; 1600];
c =[65 52 72 63 0 0 0 0]';
Beta = [5 6 7 8];
v = [1 2 3 4];

%running the file simplex which is our implementation of the simplex method
simplex(A, b, c, Beta, v);
%%
%linprog implementation
A = [1 2 2 4; 50 75 100 10; 3 1 0 0; 8 3 1 1];
c = [-65 -52 -72 -63];
b = [250 8000 4500 1600];
Aeq = [];
Beq = [];
lb = [0 0 0 0];
ub = [inf inf inf inf];
[X, Z] = linprog(c, A, b, Aeq, Beq, lb, ub);
Batches_linprog = fix(X)
cost_linprog = round(-Z)