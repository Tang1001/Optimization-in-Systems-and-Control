%% LQP Assignment 2021
%Alexandra Ministeru, Weihong Tang

clear all
%Team specific parameters
E1 = 5;
E2 = 9;
E3 = 2;

%% 1
%min f = -3.5*G - 2.2*R
f1 = [-3.5 -2.2 0 0];                                                       
Aeq1=[1 1 1 0;250 120 0 1];                                                 
beq1=[10;2000+20*E1];
lb1 = zeros(4,1);                                                           
o = optimoptions('linprog','Algorithm','dual-simplex');
[x1,fval1]=linprog(f1,[],[],Aeq1,beq1,lb1,[],o)

% Approximate integer solutions
G= ceil(x1(1));
R = floor(x1(2));
sum = G + R;
cost = 250*G + 120*R;
if (sum <= 10 && cost <= 2100) % Check if approximation fulfills constraints
    fprintf('Integer solution:\n')
    sol = [G; R]
    fprintf('Mining rate:\n')
    mining_rate = 3.5*G + 2.2*R
end

G= floor(x1(1));
R = ceil(x1(2));
sum = G + R;
cost = 250*G + 120*R;
if (sum <= 10 && cost <= 2100)
    fprintf('Integer solution:\n')
    sol = [G; R]
    fprintf('Mining rate:\n')
    mining_rate = 3.5*G + 2.2*R
end


%% 1.3
% min m = 250*G + 120*R
f13 = [250 120 0 0 0];
Aeq13=[1 1 1 0 0;-365*3.5 -365*2.2 0 1 0;250 120 0 0 1];
beq13=[10;-8000-20*E2;2000+20*E1];
lb13 = zeros(5,1);

o = optimoptions('linprog','Algorithm','dual-simplex');
[x13,fval13]=linprog(f13,[],[],Aeq13,beq13,lb13,[],o);

% Approximate and check integer solutions
G= ceil(x13(1));
R = floor(x13(2));
sum = G + R;
if (sum <= 10)
    fprintf('Integer solution:\n')
    sol = [G; R]
    fprintf('Mining rate error:\n')
    mining_error = 8180/365 - 3.5*G - 2.2*R
end

G= floor(x13(1));
R = ceil(x13(2));
sum = G + R;
if (sum <= 10)
    fprintf('Integer solution:\n')
    sol = [G; R]
    fprintf('Mining rate error:\n')
    mining_error = 8180/365 - 3.5*G - 2.2*R
end


%% 3
%x=quadprog(H,c,A,b,Aeq,beq,lb,ub,x0,options)
mm = table2array(readtable('measurements21.csv'));

q3 = [mm(1:1500,3)';mm(1:1500,4)';mm(1:1500,2)']; %3x1500
Tk3 = mm(1:1500,1);  %1500x1
Tk13 = mm(2:1501,1); %1500x1

% Hessian matrix
H = zeros(4,4);
H(1,1) = 2*Tk3'*Tk3; % 1x1 element
H(1,2:4) = 2*Tk3'*q3'; % 1x3 block
H(2:4,1) = H(1,2:4); % 3x1 block
H(2:4,2:4) = 2*q3*q3'; % 3x3 block

% Linear part
f3 = [-2*Tk13'*Tk3, -2*Tk13'*q3'];

[x3,fval3,exitflag3,output3,lambda3] = quadprog(H,f3,[],[]);

a1 = x3(4)/60;
a2 = x3(2)/60;


%% 4
%[x,fval]=linprog(f,A,b,Aeq,beq,lb,ub)
A = x3(1);
B = x3(2:4)';
T1 = 23.3;
N = 1440;

cost = mm(1:1440,5)'/1000; % phi, euro/kWh converted to euro/Wh
spent = ones(1,1440)*(0.35+0.01*E3)/1000; % psi, euro/kWh converted to euro/Wh
f4 = [(cost-spent)*(1/60),zeros(1,1440)]; % delta_t converted to hours, zeros correspond to Tk

Qout_Tamb = [mm(1:1440,4)';mm(1:1440,2)'];
Q = B(2:3)* Qout_Tamb; % B split into B(1) and B23

i = ones(N);
Aeq4 = [eye(N)*B(1),(eye(N)*(-1)+(tril(i)-tril(i,-2)-eye(N))*A)]; % matrix A: fist part corresponds to qk_in, 
                                                                  % second part to Tk
beq4 = -Q;
beq4(1) = beq4(1) - A*T1; % Add initial condition to first equation
beq4 = beq4';

lb4 = zeros(1,1440); % Lower bounds not specified for Tk
ub4 = [ones(1,1440)*125,ones(1,1440)*90]; % Constraints on maximum value of q_in and T

o = optimoptions('linprog','Algorithm','dual-simplex');
[x4,fval4,exitflag4,output4]=linprog(f4,[],[],Aeq4,beq4,lb4,ub4,[],o);

% Value of objective function at optimum
goal = cost*(x4(1:1440)+mm(1:1440,4))*(1/60)-spent*x4(1:1440)*(1/60)

% Optimized cost of mining
q_in_optimized = x4(1:1440);
q_out = mm(1:1440,4);
optimized_cost = cost*(q_out + q_in_optimized)*(1/60)