clear all;
close all;
clc;

%% Panel

c = 1; %% chord length ratio
ratio = 0.12; %% NACA0012
N = 1001; %% The number of Panels
AOA = -5*pi/180:5*pi/180:20*pi/180; %% radian

%%% Symmetric airfoil NACA0012
%%% X-Boundary points by Cosine spacing (more panels at LE & TE)
theta = linspace(0,pi,(N+1)/2+1);
for i = 1:(N+1)/2
    X_1(i) = 0.5*(1-cos(theta(i+1)));
end

theta2 = linspace(pi,0,(N+1)/2+1);
for i = (N+1)/2:N+1
    X_1(i) = 0.5*(1-cos(theta2(i-(N-1)/2)));
    if i == N+1
        X_1(i) = X_1(1);
    end
end

%%% Y-Boundary points (zero thickness at TE)
for i = 1:N+1
    Y_1(i) = 5*ratio*c*(0.2969.*sqrt(X_1(i)/c)-0.1260.*(X_1(i)/c)- ...
             0.3516.*(X_1(i)/c).^2+0.2843.*(X_1(i)/c).^3-0.1036.*(X_1(i)/c).^4);
    if i == (N+1)/2 
        Y_1(i) = 0;
    end
end

%%% NACA 2412
m = 0.02;
p = 0.4;
for i = 1:N+1
    if X_1(i) <= p*c
        theta(i) = atan((2*m/p^2)*(p-X_1(i)/c));
    else
        theta(i) = atan((2*m/(1-p)^2)*(p-X_1(i)/c));
    end
end

for i = 1:N+1
    if X_1(i) <= p*c
        Y_c(i) = (c*m/p^2)*(2*p*(X_1(i)/c)-(X_1(i)/c)^2);
    else
        Y_c(i) = (c*m/(1-p)^2)*((1-2*p)+2*p*(X_1(i)/c)-(X_1(i)/c)^2);
    end
end

for i = 1:N+1
    if i < (N+1)/2+1
        X(i) = X_1(i)-Y_1(i)*sin(theta(i));
        Y(i) = Y_c(i)+Y_1(i)*cos(theta(i));
    else 
        X(i) = X_1(i)+Y_1(i)*sin(theta(i));
        Y(i) = Y_c(i)-Y_1(i)*cos(theta(i));
    end
    if i == N+1
        X(i) = X(1);
        Y(i) = Y(1);
    end
end

%%% Control points
for i = 1:N
    x(i) = (X(i)+X(i+1))/2;
    y(i) = (Y(i)+Y(i+1))/2;
end
X_0 = [0 1];
Y_0 = [0 0];
%%% Panel figure
figure(1);
hold on
plot(X,Y,'-kx','Markersize',5);
plot(x,y,'b.','Markersize',20);
plot(X_1,Y_c,'r');
plot(X_0,Y_0,'k');
xlabel('x/c');
xlim([0 1]);
title('NACA2412 Panel N = 1001');
legend('Boundary points','Control points','Camber line','Chord line');
hold off

%% Velocity calculations

V_inf = 1;

%% Panel length

for i = 1:N
    S(i) = sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2);
end

%% Panel & Normal angles

for i = 1:N
    Phi(i) = atan2((Y(i+1)-Y(i)),(X(i+1)-X(i)));
    if Phi(i)<0
        Phi(i)=Phi(i)+2*pi;
    end
end

for i = 1:N
    Beta(i) = Phi(i)+pi/2;
    if Beta(i)>2*pi
        Beta(i) = Beta(i)-2*pi;
    end
end

%% Calculating Integral

for i=1:N
    for j=1:N
        A(i,j)=-(x(i)-X(j))*cos(Phi(j))-(y(i)-Y(j))*sin(Phi(j));
        B(i,j)=(x(i)-X(j))^2+(y(i)-Y(j))^2;
        C(i,j)=-cos(Phi(i)-Phi(j));
        D(i,j)=(y(i)-Y(j))*sin(Phi(i))+(x(i)-X(j))*cos(Phi(i));
        E(i,j)=(x(i)-X(j))*sin(Phi(j))-(y(i)-Y(j))*cos(Phi(j));
    end
end

for i=1:N
    for j=1:N
        if i==j
            J(i,j)=0;
        else
            J(i,j)=(C(i,j)/2)*log((S(j)^2+2*A(i,j)*S(j)+B(i,j))/B(i,j))+ ...
                   (D(i,j)-A(i,j)*C(i,j))/E(i,j)*(atan((S(j)+A(i,j))/E(i,j)) ...
                   -atan(A(i,j)/E(i,j)));
        end
    end
end

for j = 1:6
    for i = 1:N
        V(i,j) = sin(Phi(i)-AOA(j));
    end
end
%% Kutta condition
for j = 1:6
    J(N,:)=0;
    J(N,(N-1)/2)=1;
    J(N,(N+1)/2)=1;
    V(N,j)=0;
end
%% Vortex strength
gamma_coefficient = zeros(N,6);
gamma = zeros(N,6);
invJ = inv(J);
for j = 1:6
    gamma_coefficient(:,j) = invJ*V(:,j);
end
gamma = gamma_coefficient*(-2)*pi*V_inf;
gamma = gamma';
%% Surface Velocity, Cp calculation
V_iu = zeros(6,(N+1)/2);
Cpu = zeros(6,(N+1)/2);
for j = 1:6
    for i = 1:(N+1)/2
        V_iu(j,i) = gamma(j,i);
        Cpu(j,i) = 1-(V_iu(j,i)/V_inf)^2;
    end
end
V_il = zeros(6,(N+1)/2);
Cpl = zeros(6,(N+1)/2);
for j = 1:6
    for i = (N+1)/2:N
        V_il(j,i-(N-1)/2) = gamma(j,i);
        Cpl(j,i-(N-1)/2) = 1-(V_il(j,i-(N-1)/2)/V_inf)^2;
    end
end

figure(2);
hold on   
plot(x(1:(N+1)/2),Cpu(1,:),'r',x((N+1)/2:N),Cpl(1,:),'--r',x(1:(N+1)/2),Cpu(2,:),...
'y',x((N+1)/2:N),Cpl(2,:),'--y',x(1:(N+1)/2),Cpu(3,:),'g',x((N+1)/2:N),Cpl(3,:),'--g'...
,x(1:(N+1)/2),Cpu(4,:),'b',x((N+1)/2:N),Cpl(4,:),'--b',x(1:(N+1)/2),Cpu(5,:),...
'm',x((N+1)/2:N),Cpl(5,:),'--m',x(1:(N+1)/2),Cpu(6,:),'k',x((N+1)/2:N),Cpl(6,:),'--k');
legend('\alpha = -5, Upper','\alpha = -5, Lower','\alpha = 0, Upper', ...
'\alpha = 0, Lower','\alpha = 5, Upper','\alpha = 5, Lower',...
'\alpha = 10, Upper','\alpha = 10, Lower','\alpha = 15, Upper','\alpha = 15, Lower',...
'\alpha = 20, Upper','\alpha = 20, Lower');
grid on;
xlabel('x/c');
ylabel('Cp');
set(gca,'Ydir','reverse');
title('Cp vs x/c (1001 Panels)');
ylim([-25 2]);
xlim([0 1]);
hold off

%% Lift coefficient (from Kutta-Joukowski theorem & Cp distribution)

%Kutta-Joukowski
Circul = S*(gamma');
Cl = Circul/(1/2*V_inf*c);

%Cp distribution
for j = 1:6
    for i = 1:(N+1)/2
        if i == (N+1)/2
            Cn(j,i) = Cpl(j,i)*(X(i+(N-1)/2)-X(i+(N+1)/2));
            Ca(j,i) = Cpl(j,i)*(Y(i+(N+1)/2)-Y(i+(N-1)/2));
        else
            Cn(j,i) = Cpl(j,i)*(X(i+(N-1)/2)-X(i+(N+1)/2))-Cpu(j,i)*(X(i+1)-X(i));
            Ca(j,i) = Cpu(j,i)*(Y(i+1)-Y(i))+Cpl(j,i)*(Y(i+(N+1)/2)-Y(i+(N-1)/2));
        end
    end
Cn_sum(j) = sum(Cn(j,:))/c;
Ca_sum(j) = sum(Ca(j,:))/c;
cl(j) = Cn_sum(j)*cos(AOA(j))-Ca_sum(j)*sin(AOA(j));
end

%Cl from experiments
AOA_d = -5:5:20;
fid = fopen('Exp(NACA2412_Cl).dat','r');
data = textscan(fid, '%f%f%f', 'HeaderLines', 1, 'Collect', 1);
fclose(fid);
Cl_exp = data{1};
%Cl from vortex panel method
figure(3);
hold on
grid on
plot(Cl_exp(:,1),Cl_exp(:,2),'b',AOA_d,cl,'k');
title('Cl comparison');
legend('Experiments data','Vortex panel method','Location','SouthEast');
xlabel('Angle of attack');
ylabel('Cl');
hold off

figure(4);
hold on
grid on
x1 = -2.1464;
y1 = 0;
txt1 = '\leftarrow Cl=0 at AOA = -2.1464deg';
text(x1,y1,txt1)
plot(AOA_d,cl,'k');
title('Cl vs AOA');
legend('Vortex panel method','Location','SouthEast');
xlabel('Angle of attack');
ylabel('Cl');
hold off