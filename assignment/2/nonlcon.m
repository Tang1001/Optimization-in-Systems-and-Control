function [c,ceq] = nonlcon(vu)

ceq = [];
% Parameters
E1 = 5;
E2 = 9;
E3 = 2;

H = 5;          %tower height
g = 9.80665;
phi = 3;        %cross-flow area

%heat transfer coefficients
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
beta = 10^6;        %the penalization factor for thermal comfort

%200*5mins
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
xk_storage = zeros(4,1);

N = 144;  %number of time examples
c = zeros(N,1);

for i=0:N     %time:k=0~N
    j = i+1;    %j denotes index
    %now k=i, only use states(i)
    
    mk_dot = rho(1) * vu((N-1)+j) * phi * sqrt(2 * g * H * max(0, (xk(1)-To(j))/xk(1)));
    c(j) = -mk_dot + 0.5;
    
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
        
    elseif i == N-1 %only update the xk(1)
        xk_storage(1) = xk(1) + (qp_dot(j) + mk_dot*C(1)*(To(j) - xk(1)) ...
            + ha*A(2)*(xk(2)-xk(1)) + ha*A(3)*(xk(3)-xk(1)) ...
            + ha*A(4)*(xk(4)-xk(1))) * (dt/(rho(1)*V(1)*C(1)));
        xk(1) = xk_storage(1);
    end
end
end