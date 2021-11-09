clear all;
%% 
N = 144;  %number of time samples
t = 0:144;

%% Ex 3
% Initial conditions
vk = 0.1*ones(N-1,1);
uk = 0.1*ones(N+1,1);

lb = zeros((N-1+N+1),1);
ub = ones((N-1+N+1),1);

%%sqp
options = optimoptions('fmincon','Algorithm','sqp');
options.MaxIterations = 1000;
options.MaxFunctionEvaluations = 1.7e+06;
tic
[vu,fval,exitflag,output] = fmincon(@objective,[vk;uk],[],[],[],[],lb,ub,[],options);
time = toc;

%% Ex 4 - two different starting points
vk_diff = 0.9*ones(N-1,1);
uk_diff = 0.9*ones(N+1,1);
[vu_diff,fval_diff,exitflag_diff,output_diff] = fmincon(@objective,[vk_diff;uk_diff],[],[],[],[],lb,ub,[],options);

%plot u and v
t = 0:144;
figure(1)
plot(t(1:145),vu(144:288), 'linewidth', 1.2)
hold on
plot(t(1:145),vu_diff(144:288), 'linewidth', 1.2)
legend({'u(k)=0.1, v(k)=0.1', 'u(k)=0.9, v(k)=0.9'})
xlabel('Time index')
grid on
title('Ventilation factor u')
set(gcf,'color','w');

figure(2)
plot(t(1:143),vu(1:143), 'linewidth', 1.2)
hold on
plot(t(1:143),vu_diff(1:143), 'linewidth', 1.2)
legend({'u(k)=0.1, v(k)=0.1', 'u(k)=0.9, v(k)=0.9'})
xlabel('Time index')
grid on
title('Shading factor v')
set(gcf,'color','w');

%% plot backup energy
figure(3)
X_diff = categorical({'All starting points:0.1', 'All starting points:0.9'});
X_diff = reordercats(X_diff,{'All starting points:0.1', 'All starting points:0.9'});
y_backup = [fval fval_diff];
b= bar(X_diff,y_backup,0.4);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
ylabel('Energy (J)')
ylim([0 15e+08])
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Total backup energy, different starting points')
set(gcf,'color','w');

%% Ex 5
%%multi-start sqp
vk1 = 0*ones(N-1,1);
uk1 = 0*ones(N+1,1);
vk2 = 0.2*ones(N-1,1);
uk2 = 0.2*ones(N+1,1);
vk3 = 0.4*ones(N-1,1);
uk3 = 0.4*ones(N+1,1);
vk4 = 0.6*ones(N-1,1);
uk4 = 0.6*ones(N+1,1);
vk5 = 0.8*ones(N-1,1);
uk5 = 0.8*ones(N+1,1);
vk6 = ones(N-1,1);
uk6 = ones(N+1,1);
tic
[vu1,fval1,exitflag1,output1] = fmincon(@objective,[vk1;uk1],[],[],[],[],lb,ub,[],options);
[vu2,fval2,exitflag2,output2] = fmincon(@objective,[vk2;uk2],[],[],[],[],lb,ub,[],options);
[vu3,fval3,exitflag3,output3] = fmincon(@objective,[vk3;uk3],[],[],[],[],lb,ub,[],options);
[vu4,fval4,exitflag4,output4] = fmincon(@objective,[vk4;uk4],[],[],[],[],lb,ub,[],options);
[vu5,fval5,exitflag5,output5] = fmincon(@objective,[vk5;uk5],[],[],[],[],lb,ub,[],options);
[vu6,fval6,exitflag6,output6] = fmincon(@objective,[vk6;uk6],[],[],[],[],lb,ub,[],options);
vu_multi = (vu1+vu2+vu3+vu4+vu5+vu6)/6;

fvals = [fval1 fval2 fval3 fval4 fval5 fval6 fval];
[fval_multi,index] = min(fvals);

switch index
    case 1
        vu_multi = vu1;
    case 2
        vu_multi = vu2;
    case 3
        vu_multi = vu3;
    case 4
        vu_multi = vu4;
    case 5
        vu_multi = vu5;
    case 6
        vu_multi = vu6;
    case 7
        vu_multi = vu;
end

multi_time = toc;

%%Genetic algorithm
ga_options.MaxIterations = 5000;
ga_options.MaxFunctionEvaluations = 1.7e+06;
tic
[x_ga,fval_ga,exitflag_ga,output_ga] = ga(@objective,2*N,[],[],[],[],lb,ub,[],ga_options);
ga_time = toc;

%%Simulated annealing
an_options.MaxIterations = 50000;
an_options.MaxFunctionEvaluations = 1.7e+06;
tic
[x_sa,fval_sa,exitflag_sa,output_sa] = simulannealbnd(@objective,[vk;uk],lb,ub,an_options);
sa_time = toc;

%%Plot backup energy values
figure(4)
X_method = categorical({'Sqp algorithm','Multi-start sqp algorithm','Genetic algorithm','Simulated annealing algorithm'});
X_method = reordercats(X_method,{'Sqp algorithm','Multi-start sqp algorithm','Genetic algorithm','Simulated annealing algorithm'});
y_backup = [fval fval_multi fval_ga fval_sa];
b = bar(X_method,y_backup);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Backup energy obtained from different methods')

figure(5)
y_time = [time multi_time ga_time sa_time];
b = bar(X_method,y_time,0.4);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('Running times of different methods')

%% Ex 6 constraint mk_dot >= 0.5
[vu_ncon,fval_ncon,exitflag_ncon,output_ncon] = fmincon(@objective,[vk;uk],[],[],[],[],lb,ub,@nonlcon,options);

%Obtain the states values
%parameters
E1 = 5;
E2 = 9;
E3 = 2;

H = 5;          %tower height
g = 9.81;
phi = 3;        %cross-flow area

% heat transfer coefficients
ha = 5;         
ho = 25;
lambda = 4.5e-8;

%system properties
V = [1430; 17.5; 15; 80];
A = [0; 350; 286; 286];
C = [1005; 900; 840; 840];
rho = [1.2; 2500; 2000; 2000];
alfa = [0; 0.085; 0.2; 0.6];
tau = 0.9;
beta = 1e6;        %the penalization factor for thermal comfort

% 200*5mins
I = zeros(200+1,1);
To = zeros(200+1,1);
qp_dot = zeros(200+1,1);

I(1:47+1) = 300 + E1;
I(48+1:96+1) = 700 + E1;
I(97+1:end) = 300 + E1;

To(1:47+1) = 12 + 0.1*E2 + 273.15;
To(48+1:96+1) = 18 + 0.1*E2 + 273.15;
To(97+1:end) = 17 + 0.1*E2 + 273.15;

qp_dot(1:47+1) = 5000 + 10*E3;
qp_dot(48+1:96+1) = 25000 + 10*E3;
qp_dot(97+1:end) = 20000 + 10*E3;

dt = 5*60; %sec

Tref = 21 + 273.15; %reference temp (Kelvin)
Tsky = To - 8;      %sky temp

%initial state
Ta_0 = 16 + 273.15;
Tw1_0 = 16 + 273.15;
Tw2_0 = 16 + 273.15;
Tw3_0 = 16 + 273.15;
xk = [Ta_0; Tw1_0; Tw2_0; Tw3_0];   

N = 144;  %number of time samples

q_backup_ncon = zeros(N,1);
q_backup_ncon_sum = 0;
Ta_1 = zeros(N,1);
Tw1_1 = zeros(N-1,1);
Tw2_1 = zeros(N-1,1);
Tw3_1 = zeros(N-1,1);

for i=0:N     %time:k=0~N
    j = i+1;    %j denotes index
    %now k=i, only use states(i)
    mk_dot = rho(1) * vu_ncon((N-1)+j) * phi * sqrt(2 * g * H * max(0, (xk(1)-To(j))/xk(1)));
    q_backup = mk_dot*C(1)*abs(xk(1) - Tref)*dt + beta*(xk(1) - Tref)^2;
    
    if i<= N-2
        %update states(i+1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        
        xk_storage(2) = xk(2) + (vu_ncon(j)*I(j)*alfa(2) + lambda*(Tsky(j)^4 - xk(2)^4) ...
            + ho*(To(j) - xk(2)) + ha*(xk(2) - xk(1))) * ((A(2)*dt)/(rho(2)*V(2)*C(2)));
        
        xk_storage(3) = xk(3) + (I(j)*alfa(3) + lambda*(Tsky(j)^4 - xk(3)^4) ...
            + ho*(To(j) - xk(3)) + ha*(xk(3) - xk(1))) * ((A(3)*dt)/(rho(3)*V(3)*C(3)));
        
        xk_storage(4) = xk(4) + (0.5*vu_ncon(j)*I(j)*tau*alfa(4) ...
            + ha*(xk(4) - xk(1))) * ((A(4)*dt)/(rho(4)*V(4)*C(4)));
        
        xk = xk_storage;
        
        Ta_1(j) = xk_storage(1);
        Tw1_1(j) = xk_storage(2);
        Tw2_1(j) = xk_storage(3);
        Tw3_1(j) = xk_storage(4);
        
    elseif i == N-1 % only update xk(1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        xk(1) = xk_storage(1);
        Ta_1(j) = xk_storage(1);
    end
        
    q_backup_ncon(j) = q_backup;
    q_backup_ncon_sum = q_backup_ncon_sum + q_backup;
end

% plot uk and vk
figure(6)
plot(t(1:145),vu_ncon(144:288), 'linewidth', 1.2)
hold on
plot(t(1:143),vu_ncon(1:143), 'linewidth', 1.2)
legend({'Ventilation factor u', 'Shading factor v'})
xlabel('Time index')
grid on
title('Ventilation factor u and Shading factor v with the air flow rates constraint')
set(gcf,'color','w');

% plot states
t = 0:144;
figure(7)
subplot(221)
plot(t,[Ta_0; Ta_1]-273.15,'linewidth',1.1)
hold on
plot(t,To(1:145)-273.15,'linewidth',1.1)
hold on
plot(t,Tref*ones(145,1)-273.15,'linewidth',1.1)
xlabel('Time index')
ylabel('Temperature (degC)')
legend({'T_{a}','T_{o}','T_{ref}'},'location','southeast')
grid on
title('T_{a}')
set(gcf,'color','w');

subplot(222)
plot(t(1:144),[Tw1_0; Tw1_1]-273.15,'linewidth',1.1)
xlabel('Time index')
ylabel('Temperature (degC)')
grid on
title('T_{w1}')
set(gcf,'color','w');

subplot(223)
plot(t(1:144),[Tw2_0; Tw2_1]-273.15,'linewidth',1.1)
xlabel('Time index')
ylabel('Temperature (degC)')
grid on
title('T_{w2}')
set(gcf,'color','w');

subplot(224)
plot(t(1:144),[Tw3_0; Tw3_1]-273.15,'linewidth',1.1)
xlabel('Time index')
ylabel('Temperature (degC)')
grid on
title('T_{w3}')
set(gcf,'color','w');

% m_dot1
figure(8)
plot(t,q_backup_ncon, 'linewidth', 1.2)
xlabel('Time index')
ylabel('Mass flow (kg/s)')
grid on
title('Mass flow m with the air flow rates constraint')
set(gcf,'color','w');


%% Ex 7 
%plot uk and vk
figure(9)
plot(t(1:145),vu(144:288),'linewidth', 1.2)
hold on
plot(t(1:145),vu_diff(144:288), 'linewidth', 1.2)
hold on
plot(t(1:145),vu_ncon(144:288), 'linewidth', 1.2)
hold on
plot(t(1:145),vu_multi(144:288), 'linewidth', 1.2)
hold on
plot(t(1:145),x_ga(144:288), 'linewidth', 1.2)
hold on
plot(t(1:145),x_sa(144:288), 'linewidth', 1.2)
hold on
legend({'Sqp algorithm', 'Sqp algorithm(using different starting points)','Adding the air flow rates constraints','Multi-start sqp algorithm','Genetic algorithm','Simulated annealing algorithm'})
xlabel('Time index')
grid on
title('Ventilation factor u obtained in Task 2-6')
set(gcf,'color','w');

figure(10)
plot(t(1:143),vu(1:143),'linewidth', 1.2)
hold on
plot(t(1:143),vu_diff(1:143), 'linewidth', 1.2)
hold on
plot(t(1:143),vu_ncon(1:143), 'linewidth', 1.2)
hold on
plot(t(1:143),vu_multi(1:143), 'linewidth', 1.2)
hold on
plot(t(1:143),x_ga(1:143), 'linewidth', 1.2)
hold on
plot(t(1:143),x_sa(1:143), 'linewidth', 1.2)
hold on
legend({'Sqp algorithm', 'Sqp algorithm(using different starting points)','Adding the air flow rates constraints','Multi-start sqp algorithm','Genetic algorithm','Simulated annealing algorithm'})
xlabel('Time index')
grid on
title('Shading factor v obtained in Task 2-6')
set(gcf,'color','w');

%% q_backup

q_backup1 = zeros(N,1);
q_backup1_sum = 0;
xk = [Ta_0; Tw1_0; Tw2_0; Tw3_0];

for i=0:N     %time:k=0~N
    j = i+1;    %j denotes index
    %now k=i, only use states(i)
    mk_dot = rho(1) * vu((N-1)+j) * phi * sqrt(2 * g * H * max(0, (xk(1)-To(j))/xk(1)));
    q_backup = mk_dot*C(1)*abs(xk(1) - Tref)*dt + beta*(xk(1) - Tref)^2;
    
    if i<= N-2
        %update states(i+1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        
        xk_storage(2) = xk(2) + (vu(j)*I(j)*alfa(2) + lambda*(Tsky(j)^4 - xk(2)^4) ...
            + ho*(To(j) - xk(2)) + ha*(xk(2) - xk(1))) * ((A(2)*dt)/(rho(2)*V(2)*C(2)));
        
        xk_storage(3) = xk(3) + (I(j)*alfa(3) + lambda*(Tsky(j)^4 - xk(3)^4) ...
            + ho*(To(j) - xk(3)) + ha*(xk(3) - xk(1))) * ((A(3)*dt)/(rho(3)*V(3)*C(3)));
        
        xk_storage(4) = xk(4) + (0.5*vu(j)*I(j)*tau*alfa(4) ...
            + ha*(xk(4) - xk(1))) * ((A(4)*dt)/(rho(4)*V(4)*C(4)));
        
        xk = xk_storage;
        
    elseif i == N-1 % only update xk(1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        xk(1) = xk_storage(1);
        Ta_1(j) = xk_storage(1);
    end
        
    q_backup1(j) = q_backup;
    q_backup1_sum = q_backup1_sum + q_backup;
end

q_backup_nocontrol = zeros(N,1);
q_backup_nocontrol_sum = 0;
xk = [Ta_0; Tw1_0; Tw2_0; Tw3_0];
vu = ones(2*N);

for i=0:N     %time:k=0~N
    j = i+1;    %j denotes index
    %now k=i, only use states(i)
    mk_dot = rho(1) * vu((N-1)+j) * phi * sqrt(2 * g * H * max(0, (xk(1)-To(j))/xk(1)));
    q_backup = mk_dot*C(1)*abs(xk(1) - Tref)*dt + beta*(xk(1) - Tref)^2;
    
    if i<= N-2
        %update states(i+1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        
        xk_storage(2) = xk(2) + (vu(j)*I(j)*alfa(2) + lambda*(Tsky(j)^4 - xk(2)^4) ...
            + ho*(To(j) - xk(2)) + ha*(xk(2) - xk(1))) * ((A(2)*dt)/(rho(2)*V(2)*C(2)));
        
        xk_storage(3) = xk(3) + (I(j)*alfa(3) + lambda*(Tsky(j)^4 - xk(3)^4) ...
            + ho*(To(j) - xk(3)) + ha*(xk(3) - xk(1))) * ((A(3)*dt)/(rho(3)*V(3)*C(3)));
        
        xk_storage(4) = xk(4) + (0.5*vu(j)*I(j)*tau*alfa(4) ...
            + ha*(xk(4) - xk(1))) * ((A(4)*dt)/(rho(4)*V(4)*C(4)));
        
        xk = xk_storage;
        
    elseif i == N-1 % only update xk(1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        xk(1) = xk_storage(1);
    end
        
    q_backup_nocontrol(j) = q_backup;
    q_backup_nocontrol_sum = q_backup_nocontrol_sum + q_backup;
end

figure(11)
plot(t,q_backup1, 'linewidth', 1.2)
hold on
plot(t,q_backup_ncon, 'linewidth', 1.2)
hold on
plot(t,q_backup_nocontrol, 'linewidth', 1.2)
legend({'Unconstrained air flow', 'Constrained air flow','No control'})
xlabel('Time index')
ylabel('Energy (J)')
grid on
title('Demanded energy q_{backup}')
set(gcf,'color','w');

