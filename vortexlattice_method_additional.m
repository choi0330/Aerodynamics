clear all;
close all;
clc;

%% Setting up points

% 1. (b)
AR = 5;
S = 1; % 날개 넓이 1이라 가정
b = sqrt(S*AR);
alpha = 10; %deg
V_inf = 1; % 자유류 속도 1이라 가정

% Bound points

for i = 1:11
    x(i) = 1/AR*b/4;
    y(i) = -0.5*b+0.1*b*(i-1);
end


% Control points
for i = 1:10
    xc(i) = 1/AR*b*3/4;
    yc(i) = -0.45*b+0.1*b*(i-1);
end

%% Calculating downwash & normal velocity
for i = 1:10
    for j = 1:10
        W(i,j) = 1/(4*pi)*((1/((xc(i)-x(j))*(yc(i)-y(j+1))-(xc(i)-x(j+1))*(yc(i)-...
        y(j))))*(((x(j+1)-x(j))*(xc(i)-x(j))+(y(j+1)-y(j))*(yc(i)-y(j)))/sqrt((xc(i)-...
        x(j))^2+(yc(i)-y(j))^2)-((x(j+1)-x(j))*(xc(i)-x(j+1))+(y(j+1)-y(j))*(yc(i)-...
        y(j+1)))/sqrt((xc(i)-x(j+1))^2+(yc(i)-y(j+1))^2))+(1/(y(j)-yc(i)))*(1+...
        (xc(i)-x(j))/sqrt((xc(i)-x(j))^2+(yc(i)-y(j))^2))-(1/(y(j+1)-yc(i)))*(1+...
        (xc(i)-x(j+1))/sqrt((xc(i)-x(j+1))^2+(yc(i)-y(j+1))^2)));
    end
end
        
for i = 1:10
    V(i)=-V_inf*sin(alpha*pi/180); %rad
end

gamma = inv(W)*V';
Cl = 2*sum(gamma)*0.1*b/(V_inf*S);
