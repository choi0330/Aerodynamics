clear all;
close all;
clc;

%% Setting up points

% 1. (a)
AR = 10;
S = 1; % 날개 넓이 1이라 가정
b = sqrt(S*AR);
alpha = 5; %deg
V_inf = 1; % 자유류 속도 1이라 가정
Lambda = 0:1:30;

% Bound points
for j = 1:31
    for i = 1:11
        if i < 6
            x(i,j) = 1/AR*b/4+(0.6-i/10)*b*tan(Lambda(j)*pi/180);
            y(i,j) = -0.5*b+0.1*b*(i-1);
        elseif i == 6
            x(i,j) = 1/AR*b/4;
            y(i,j) = -0.5*b+0.1*b*(i-1);
        elseif i > 6
            x(i,j) = 1/AR*b/4+(i/10-0.6)*b*tan(Lambda(j)*pi/180);
            y(i,j) = -0.5*b+0.1*b*(i-1);
        end
    end
end

% Control points
for j = 1:31
    for i = 1:10
        if i < 6
            xc(i,j) = 1/AR*b*3/4+(0.55-i/10)*b*tan(Lambda(j)*pi/180);
            yc(i,j) = -0.45*b+0.1*b*(i-1);
        elseif i >= 6
            xc(i,j) = 1/AR*b*3/4+(-0.55+i/10)*b*tan(Lambda(j)*pi/180);
            yc(i,j) = -0.45*b+0.1*b*(i-1);
        end
    end
end

%% Calculating downwash & normal velocity
for k = 1:31
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

gamma = zeros(10,31);
for k = 1:31
    gamma(:,k) = inv(W(:,:,k))*V';
end

for k = 1:31
    Cl(k) = 2*sum(gamma(:,k))*0.1*b/(V_inf*S);
end

%% Plotting
figure(1);
plot(Lambda,Cl);
grid on;
xlabel('Swept-back angle(deg)');
ylabel('Cl');
title('Cl vs Swept-back angle (\alpha=5deg, AR=10)');
legend('Vortex lattice method(VLM)');
