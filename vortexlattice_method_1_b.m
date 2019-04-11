clear all;
close all;
clc;

%% Setting up points

% 1. (b)
AR = 5:0.5:10;
S = 1; % 날개 넓이 1이라 가정
b = sqrt(S*AR);
alpha = 5; %deg
V_inf = 1; % 자유류 속도 1이라 가정

% Bound points
for j = 1:11
    for i = 1:11
        x(i,j) = 1/AR(j)*b(j)/4;
        y(i,j) = -0.5*b(j)+0.1*b(j)*(i-1);
    end
end

% Control points
for j = 1:11
    for i = 1:10
        xc(i,j) = 1/AR(j)*b(j)*3/4;
        yc(i,j) = -0.45*b(j)+0.1*b(j)*(i-1);
    end
end

%% Calculating downwash & normal velocity
for k = 1:11
    for i = 1:10
        for j = 1:10
        W(i,j,k) = 1/(4*pi)*((1/((xc(i,k)-x(j,k))*(yc(i,k)-y(j+1,k))-(xc(i,k)-x(j+1,k))*(yc(i,k)-...
        y(j,k))))*(((x(j+1,k)-x(j,k))*(xc(i,k)-x(j,k))+(y(j+1,k)-y(j,k))*(yc(i,k)-y(j,k)))/sqrt((xc(i,k)-...
        x(j,k))^2+(yc(i,k)-y(j,k))^2)-((x(j+1,k)-x(j,k))*(xc(i,k)-x(j+1,k))+(y(j+1,k)-y(j,k))*(yc(i,k)-...
        y(j+1,k)))/sqrt((xc(i,k)-x(j+1,k))^2+(yc(i,k)-y(j+1,k))^2))+(1/(y(j,k)-yc(i,k)))*(1+...
        (xc(i,k)-x(j,k))/sqrt((xc(i,k)-x(j,k))^2+(yc(i,k)-y(j,k))^2))-(1/(y(j+1,k)-yc(i,k)))*(1+...
        (xc(i,k)-x(j+1,k))/sqrt((xc(i,k)-x(j+1,k))^2+(yc(i,k)-y(j+1,k))^2)));
        end
    end
end
        
for i = 1:10
        V(i)=-V_inf*sin(alpha*pi/180); %rad
end

gamma = zeros(10,11);
for k = 1:11
    gamma(:,k) = inv(W(:,:,k))*V';
end

for k = 1:11
    Cl(k) = 2*sum(gamma(:,k))*0.1*b(k)/(V_inf*S);
end

%% Plotting
figure(1);
plot(AR,Cl,'Linewidth',2);
grid on;
xlabel('AR');
ylabel('Cl');
title('Cl vs Aspect Ratio (\alpha=5deg)');
legend('Vortex lattice method(VLM)');
