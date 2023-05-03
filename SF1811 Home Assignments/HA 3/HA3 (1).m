%% SF1811, HA 3 Jakob Amaya Scott, Alexander Råberg
%To run the different plots, you might need to remove the "%" signs.


%exercise 1

r = [];

%creating r from 3 to 9 with step length 0.25
for i = 3:0.25:9
    r(end+1) = i;
end

%given code from assignment
n=8;
Corr=zeros(n,n);
for i=1:n
    for j=1:n
        Corr(i,j)=(-1)^abs(i-j)/(abs(i-j)+1);
    end
end
sigma=zeros(n,1);
mu=zeros(n,1);
sigma(1)=2;
mu(1)=3;
for i=1:n-1
    sigma(i+1)=sigma(i)+2*rand;
    mu(i+1)=mu(i)+1;
end
D=diag(sigma);
C2=D*Corr*D;
C=0.5*(C2+C2');

%setting up the problem to fit the quadprog in matlab
f = [];
A = [];
H = 2*C;
b = [];
%creating a column vector with one with the same size as mu, (8 columns, 1 row)
I = ones(1,size(mu,1));
%equality constraint Aeq
Aeq = [mu';I];
%lower bound lb, with zeros
lb = zeros(1,n);
%upper bound 
ub = [inf inf inf inf inf inf inf inf];
new_sigma = [];
i = 1;

while i <= size(r,2)
    %equality constraint Beq
    Beq = [r(i),1];
    [X,fV] = quadprog(H,f,A,b,Aeq,Beq,lb,ub);
    sigma = sqrt(fV);
    new_sigma(end+1) = sigma;
    i = i + 1;
end
hold on
xlabel('Sigma vector values')
ylabel('Mu, Expected return')
title('Excerice 1 & 2')
plot(new_sigma,r,'o', 'LineWidth', 4, 'Color', 'Blue')


%excercise 2
new_sigma_2 = [];
A = [I];
b = 1;
f = [];
H = 2*C;
I = ones(1,size(mu,1));
Aeq = [mu'];
lb = zeros(1,n);
ub = [inf inf inf inf inf inf inf inf];
new_sigma_2 = []';
i = 1;
while i <= size(r,2)
    Beq = [r(i)];
    [X,fV] = quadprog(H,f,A,b,Aeq,Beq,lb,ub);
    sigma = sqrt(fV);
    new_sigma_2(end+1) = sigma;
    i = i + 1;
end
% hold on
plot(new_sigma_2,r,'*', 'LineWidth', 4, 'Color', 'red')
% legend({'y = Sigma vector 1','y = Sigma vector with possibility to save the not invested fraction'},'Location','northwest')
hold off


% excericise 3 
%mu'*x >= r =>
%mu'*x = -r


new_sigma_3 = [];
A = [-mu'];
Beq = 1;
f = [];
H = 2*C;
I = ones(1,size(mu,1));
Aeq = [I];
lb = zeros(1,n);
ub = [inf inf inf inf inf inf inf inf];
new_sigma_3 = []';
new_mu = [];

i = 1;

while i <= size(r,2)
    b = [-r(i)];
    [X,fV] = quadprog(H,f,A,b,Aeq,Beq,lb,ub);
    new_mu(end+1) = mu'*X;
    sigma = sqrt(fV);
    new_sigma_3(end+1) = sigma;
    i = i + 1;
end

%hold on
%plot(new_sigma_3,new_mu,'*', 'LineWidth', 4, 'Color', 'red')
%legend({'y = Sigma vector 1','y = Sigma vector with possibility to save the not invested fraction'},'Location','northwest')
%hold off


%excericesie 4
new_sigma_4 = [];
A = [];
b = [];
f = [];
H = 2*C;
I = ones(1,size(mu,1));
Aeq = [mu';I];
lb = zeros(1,n);
ub = [inf inf inf inf inf inf inf inf];
new_mu_2 = [];
i = 1;

while i <= size(r,2)
    Beq = [r(i);1];
    [X,fV] = quadprog(H,f,A,b,Aeq,Beq);
    new_mu_2(end+1) = mu'*X;
    sigma = sqrt(fV);
    new_sigma_4(end+1) = sigma;
    i = i + 1;
end

% hold on
% plot(new_sigma_4,new_mu_2,'*', 'LineWidth', 4, 'Color', 'red')
% hold off








