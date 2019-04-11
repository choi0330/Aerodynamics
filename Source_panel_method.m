clear all;
close all;
clc;

%% Panel

c = 1; %% chord length ratio
ratio = 0.12; %% NACA0012
N = 100; %% The number of Panels


%%% X-Boundary points
for i = 1:N/2
    X(i) = (i-1)/(N/2)*c;
end
for i = N/2+1:N+1
    X(i) = c-(i-N/2-1)/(N/2)*c;
end

%%% Cosine spacing (more panels at LE)
theta = linspace(0,pi,N/2+1);
for i = 1:N/2+1
    X(i) = 0.5*X(i)*(1-cos(theta(i)));
end
theta2 = linspace(pi,0,N/2+1);
for i = N/2+1:N+1
    X(i) = 0.5*X(i)*(1-cos(theta2(i-N/2)));
end

%%% Y-Boundary points (zero thickness at TE)
for i = 1:N+1
    Y(i) = 5*ratio*c*(0.2969.*sqrt(X(i)/c)-0.1260.*(X(i)/c) ...
           -0.3516.*(X(i)/c).^2+0.2843.*(X(i)/c).^3-0.1036.*(X(i)/c).^4);
    if i > N/2+1
        Y(i) = -5*ratio*c*(0.2969.*sqrt(X(i)/c)-0.1260.*(X(i)/c) ...
               -0.3516.*(X(i)/c).^2+0.2843.*(X(i)/c).^3-0.1036.*(X(i)/c).^4);
    end
end

%%% Control points
for i = 1:N
    x(i) = (X(i)+X(i+1))/2;
    y(i) = (Y(i)+Y(i+1))/2;
end

%%% Panel figure
figure(1);
hold on
plot(X,Y,'-kx','Markersize',5);
plot(x,y,'b.','Markersize',20);
xlabel('x/c');
title('NACA0012 Panel N = 100');
legend('Boundary points','Control points');
hold off

%% Velocity calculations

V_inf = 1;

%%% Panel length

for i = 1:N
    S(i) = sqrt((X(i+1)-X(i)).^2+(Y(i+1)-Y(i))^2);
end

%%% Panel & Normal angles

for i = 1:N
    Phi(i) = atan2((Y(i+1)-Y(i)),(X(i+1)-X(i)));
    if Phi(i)<0
        Phi(i)=Phi(i)+2*pi;
    end
end

for i = 1:N
    Beta(i) = Phi(i)+pi/2;
end

%%% Calculating Integral

for i = 1:N
    for j = 1:N
        A(i,j) = -(x(i)-X(j))*cos(Phi(j))-(y(i)-Y(j))*sin(Phi(j));
        B(i,j) = (x(i)-X(j))^2+(y(i)-Y(j))^2;
        C(i,j) = sin(Phi(i)-Phi(j));
        D(i,j) = (y(i)-Y(j))*cos(Phi(i))-(x(i)-X(j))*sin(Phi(i));
        E(i,j) = (x(i)-X(j))*sin(Phi(j))-(y(i)-Y(j))*cos(Phi(j));
    end
end

for i = 1:N
    for j = 1:N
        if i==j
            I(i,j) = pi;
        else
            S(j) = sqrt((X(j+1)-X(j))^2+(Y(j+1)-Y(j))^2);
            I(i,j) = (C(i,j)/2)*log((S(j)^2+2*A(i,j)*S(j)+B(i,j))/B(i,j)) ...
            +((D(i,j)-A(i,j)*C(i,j))/E(i,j))*(atan((S(j)+A(i,j))/E(i,j)) ...
            -atan(A(i,j)/E(i,j)));
        end
    end
end

for i = 1:N
    V(i) = sin(Phi(i));
end
V=V';

%%% Source strength
lambda_coefficient = inv(I)*V;
lambda = lambda_coefficient*2*pi*V_inf;
sum = 0;

figure(2);
hold on
plot(x(1:50),lambda(1:50),'ko',x(51:100),lambda(51:100),'k*');
grid on;
xlim([-0.1 1]);
xlabel('Control points - x/c');
ylabel('lambda i');
legend('Upper surface source strength','Lower surface source strength');
title('Source strength lambda');
hold off

for i = 1:N
    accu(i) = (S(i)*lambda(i));
    sum = accu(i)+sum;
end

accuracy = sum;

%%% Surface Velocity

for i = 1:N
    for j = 1:N
        if i == j
            J(i,j) = 0;
        else
            J(i,j) = (D(i,j)-A(i,j)*C(i,j))/(2*E(i,j))*log((S(j)^2 ...
            +2*A(i,j)*S(j)+B(i,j))/B(i,j))-C(i,j)*(atan((S(j) ...
            +A(i,j))/E(i,j))-atan(A(i,j)/E(i,j)));
        end
    end
end

V_s=J*lambda/(2*pi);

%% Cp calculation

for i = 1:N
    V_i(i) = V_inf*cos(Phi(i))+V_s(i);
    Cp(i) = 1-(V_i(i)/V_inf)^2;
end

figure(3);
hold on
plot(x(1:50),Cp(1:50),'k');
plot(x(51:100),Cp(51:100),'k*');
legend('Upper surface','Lower surface');
grid on;
xlabel('x/c');
ylabel('Cp');
xlim([-0.1 1.1]);
ylim([-0.5 1.1]);
set(gca,'XTick',-0.1:0.1:1.1);
set(gca,'YTick',-0.5:0.1:1.1);
set(gca,'Ydir','reverse');
title('Cp vs x/c');
hold off

%% Flow around the airfoil

%%% infinitesimally vertical small panel (normal vector direction : +x)

xp = (-0.5:0.05:1.5);
yp = (-0.35:0.01:0.35);
Phi_p = 3*pi/2;

figure(4);
hold on;
for jj = 1:length(yp)
    for ii = 1:length(xp)
        for j = 1:N
            A_p(j) = -(xp(ii)-X(j))*cos(Phi(j))-(yp(jj)-Y(j))*sin(Phi(j));
            B_p(j) = (xp(ii)-X(j))^2+(yp(jj)-Y(j))^2;
            C_p(j) = sin(Phi_p-Phi(j));
            D_p(j) = (yp(jj)-Y(j))*cos(Phi_p)-(xp(ii)-X(j))*sin(Phi_p);
            E_p(j) = (xp(ii)-X(j))*sin(Phi(j))-(yp(jj)-Y(j))*cos(Phi(j));
            I_p(j) = (C_p(j)/2)*log((S(j)^2+2*A_p(j)*S(j)+B_p(j))/B_p(j)) ...
            +((D_p(j)-A_p(j)*C_p(j))/E_p(j))*(atan((S(j)+A_p(j))/E_p(j)) ...
            -atan(A_p(j)/E_p(j)));
            J_p(j) = (D_p(j)-A_p(j)*C_p(j))/(2*E_p(j))*log((S(j)^2 ...
            +2*A_p(j)*S(j)+B_p(j))/B_p(j))-C_p(j)*(atan((S(j) ...
            +A_p(j))/E_p(j))-atan(A_p(j)/E_p(j)));
            V_xp(j) = I_p(j)*lambda(j)/(2*pi);
            V_yp(j) = -J_p(j)*lambda(j)/(2*pi);
        end
        %%% outside the airfoil
        V_x = 0;
        for i = 1:N
            V_x = V_x+V_xp(i);
        end
        V_x = V_x+V_inf;
        V_y = 0;
        for i = 1:N
            V_y = V_y+V_yp(i);
        end
        %%% inside the airfoil
        for i = 1:N/2+1
            y_temp(i) = 5*ratio*c*(0.2969.*sqrt(xp(ii)/c)-0.1260.*(xp(ii)/c) ...
            -0.3516.*(xp(ii)/c).^2+0.2843.*(xp(ii)/c).^3-0.1036.*(xp(ii)/c).^4);
            if xp(ii) > 0 && xp(ii) < c
                if yp(jj) < y_temp(i) && yp(jj) >= 0
                    V_x = 0;
                    V_y = 0;
                end
                if yp(jj) > -y_temp(i) && yp(jj) <= 0
                    V_x = 0;
                    V_y = 0;
                end
            end
        end
        quiver(xp(ii),yp(jj),V_x,V_y,0.05,'k');
    end
end
plot(X,Y,'b');
xlim([-0.5 1.5]);
ylim([-0.35 0.35]);
legend('Velocity Vector');
title('Velocity field around the airfoil');
hold off
