%% Cilindricas
N = 100;

s = linspace(0,50,N)
theta = linspace(-pi,pi,N)
E = 1.5;
R = 2; 
for i= 1:length(s)
    for j = 1:length(theta)
        ro(i,j) = (((E * (R^2))/s(i)) - E*s(i))*cos(theta(j));
    end 
end
[x y z] = pol2cart(theta,s,ro);
[X Y] = meshgrid(x,y);
figure(1)
surface(X,Y,z)
title('Cilindricas')
zlabel('\rho(s,\theta)');
view(3)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Esféricas caso 1
N = 100;
r = linspace(0,50,N);
theta = linspace(0,pi,N);
E = 1.5;
R = 2; 

for i= 1:length(r)
    for j = 1:length(theta)
        V(i,j) =   (((-663 * R^5)/(280*E))/(r(i)^4))*(0.5*(5*(cos(theta(j)))^3 - 3*cos(theta(j)))); 
    end
end
[x y z] = sph2cart(theta,r,V);
[X Y] = meshgrid(x,y);
%Graficar
figure(2)
surface(x,y,z)
title('Esfericas Caso 1');
xlabel('x');ylabel('y');
zlabel('V(r,\theta)');

view(3)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Esféricas caso 2 
N = 100;
r = linspace(0,50,N);
theta = linspace(0,pi,N);
E = 1.5;
R = 2; 
k = 1.6;
Vo = 5;

for i= 1:length(r)
    for j = 1:length(theta)
        if r(i) < R
            V(i,j) =   Vo + ((k*8*(r(i))^2)/(E*R*15))* ((3/2)*(cos(theta(j)))^2 - 1/2); 
        end
        if r(i) >= R
            V(i,j) =   (k*(R^4)*8/(E*15*(r(i))^3))* ((3/2)*(cos(theta(j)))^2 - 1/2); 
        end
    end
end
[x y z] = sph2cart(theta,r,V);
[X Y] = meshgrid(x,y);
%Graficar
figure(3)
surface(x,y,z)
title('Esfericas Caso 2');
xlabel('x');ylabel('y');
zlabel('V(r,\theta)');

view(3)
hold off

%% Campo Eléctrico esfericas caso 1
N = 100;
r = linspace(0,50,N);
theta = linspace(0,pi,N);
E = 1.5;
R = 10; 

for i= 1:length(r)
    for j = 1:length(theta)
        if r(i) < R
            V(i,j) = ((-663*3*r(i)^2)/(280*E*R^2))*(0.5*(5*(cos(theta(j)))^3 - 3*cos(theta(j))));
        end
        if r(i) >= R
            V(i,j) = ((663 *4 * R^5)/(280*E*r(i)^5))*(0.5*(5*(cos(theta(j)))^3 - 3*cos(theta(j))));
        end
    end
end
[x y z] = sph2cart(theta,r,V);
[X Y] = meshgrid(x,y);
%Graficar
figure(4)
surface(x,y,z)
title('Campo eléctrico esfericas caso 1');
xlabel('x');ylabel('y');
zlabel('E(r,\theta)');

view(3)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Campo Eléctrico esfericas caso 2 
N = 100;
r = linspace(0,50,N);
theta = linspace(0,pi,N);
E = 1.5;
R = 10; 
k = 1.6;
Vo = 5;

for i= 1:length(r)
    for j = 1:length(theta)
        if r(i) < R
            V(i,j) =   ((2*k*8*r(i))/(E*R*15))* ((3/2)*(cos(theta(j)))^2 - 1/2); 
        end
        if r(i) >= R
            V(i,j) =   -Vo*R/r(i)^2 - (3*k*(R^4)*8/(E*15*(r(i))^4))* ((3/2)*(cos(theta(j)))^2 - 1/2); 
        end
    end
end
[x y z] = sph2cart(theta,r,V);
%Graficar
figure(5)
surface(x,y,z)
title('Campo eléctrico esféricas caso 2');
xlabel('x');ylabel('y');
zlabel('E(r,\theta)');

view(3)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cartesianas caso 1.a 
N = 100;
x = linspace(0,10,N);
y = linspace(0,10,N);
a = 5;
b = 5;
n = 2; % 2,5,10,20
w = n*pi()/b;

fun = @(l) atan(l/a).*sin(w*l);
q = integral(fun,0,b);

for i= 1:length(x)
    for j = 1:length(y)
        V(i,j) =   (2*sinh(x(i)*w)*sinh(y(j)*w)*q)/(b*sinh(w*a)); 
    end
end

[X Y] = meshgrid(x,y);
%Graficar
figure(6)
surface(x,y,V)
xlabel('x');ylabel('y');zlabel('V(x,y)');
title('Cartesianas caso 1.a')
view(3)
hold off
%% Cartesianas caso 1.b 
N = 100;
x = linspace(0,10,N);
y = linspace(0,10,N);
a = 5;
b = 5;
n = 2; % 2,5,10,20
w = n*pi()/a;
fun = @(l) sin(l*w).*(2*(l.^3) + 5);
        q = integral(fun,0,a);

for i= 1:length(x)
    for j = 1:length(y) 
        V(i,j)=(2/a)*q*((cosh(w*x(i))*sin(w*y(j))/cosh(w*b)));
    end
end

[X Y] = meshgrid(x,y);
%Graficar
figure(7)
surface(x,y,V)
xlabel('x');ylabel('y');zlabel('V(x,y)');
title('Cartesianas caso 1.b')
view(3)
hold off

%% Cartesianas caso 2 (Convergencia de sumatorias)
N = 100;
x = linspace(0,10,N);
y = linspace(0,10,N);
a = 5;
b = 5;

for i= 1:length(x)
    for j = 1:length(y)
        
        
        
        
        
        
        
        
        
        [F p] = sumatoria(x(i),y(i),a);
        V(i,j) = F; 
    end
end

[X Y] = meshgrid(x,y);
% Graficar
figure(8)
surface(X,Y,V)
xlabel('x');ylabel('y');zlabel('V(x,y)');
title('Cartesianas caso 2 (Convergencia de la serie')
view(3)
hold off

%% Cartesianas caso 2 
N = 100;
x = linspace(0,10,N);
y = linspace(0,10,N);
a = 5;
b = 5;
n = 2; % 2,5,10,20
w = n*pi()/a;
%T
T = 2*a/n;
%Fn
fun = @(l) (sin(l*w)).^2;
q = integral(fun,0,T);
Fn = (2/T)*q;

for i= 1:length(x)
    for j = 1:length(y) 
        V(i,j) =  Fn*exp(w*x(i))*sin(w*y(i)); 
    end
end

[X Y] = meshgrid(x,y);
% Graficar
figure(9)
surface(X,Y,V)
xlabel('x');ylabel('y');zlabel('V(x,y)');
title('Cartesianas caso 2')
view(3)
hold off

%% Cartesianas caso 3 
N = 100;
x = linspace(0,10,N);
y = linspace(0,10,N);
a = 5;
b = 5;
n = 2; % 2,5,10,20
m = 2;
w = n*pi()/a;
f = m*pi/b;


fun = @(l) sin(w*l).*(l.^2)*(4/(a*b));
funa = @(r) sin(f*r).*r*1;
q = integral(fun,0,b)
u = integral(funa,0,a)

for i= 1:length(x)
    for j = 1:length(y) 
        V(i,j) =  q*u*sin(w*y(j))*sin(f*x(i)); 
    end
end

[X Y] = meshgrid(x,y);
% Graficar
figure(10)
surface(X,Y,V)
xlabel('x');ylabel('y');zlabel('V(x,y)');
title('Cartesianas caso 3')
view(3)
hold off
