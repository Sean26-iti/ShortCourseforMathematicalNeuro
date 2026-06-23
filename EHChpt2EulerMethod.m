T = 5;
delta_t = 0.0001;
omega =10;
N=T/delta_t;
x1 = zeros(1,N+1);
x2 = zeros(1,N+1);
x1(1)=0.1;
x2(1)=0;
x1_dot=@(x1,x2,omega)x1-omega*x2-(x1.^2+x2^2)*x1;
x2_dot=@(x1,x2,omega)omega*x1+x2-(x1.^2+x2^2)*x2;
for i =1:N
    x1(i+1)=x1(i) + delta_t * x1_dot(x1(i),x2(i),omega);
    x2(i+1)=x2(i) + delta_t * x2_dot(x1(i),x2(i),omega);
end
figure;
plot(x1,x2)
xlabel('x_1')
ylabel('x_2')
title('Phase Plane Plot, step=0.0001')
grid on
axis equal
delta_t=0.001;
N=T/delta_t;
x1 = zeros(1,N+1);
x2 = zeros(1,N+1);
x1(1)=0.1;
x2(1)=0;
for i =1:N
    x1(i+1)=x1(i) + delta_t * x1_dot(x1(i),x2(i),omega);
    x2(i+1)=x2(i) + delta_t * x2_dot(x1(i),x2(i),omega);
end
figure;
plot(x1,x2)
xlabel('x_1')
ylabel('x_2')
title('Phase Plane Plot, step=0.001')
grid on
axis equal
delta_t=0.01;
N=T/delta_t;
x1 = zeros(1,N+1);
x2 = zeros(1,N+1);
x1(1)=0.1;
x2(1)=0;
for i =1:N
    x1(i+1)=x1(i) + delta_t * x1_dot(x1(i),x2(i),omega);
    x2(i+1)=x2(i) + delta_t * x2_dot(x1(i),x2(i),omega);
end
figure;
plot(x1,x2)
xlabel('x_1')
ylabel('x_2')
title('Phase Plane Plot, step=0.01')
grid on
axis equal
delta_t=0.1;
N=T/delta_t;
x1 = zeros(1,N+1);
x2 = zeros(1,N+1);
x1(1)=0.1;
x2(1)=0;
for i =1:N
    x1(i+1)=x1(i) + delta_t * x1_dot(x1(i),x2(i),omega);
    x2(i+1)=x2(i) + delta_t * x2_dot(x1(i),x2(i),omega);
end
figure;
plot(x1,x2)
xlabel('x_1')
ylabel('x_2')
title('Phase Plane Plot, step=0.1')
grid on
axis equal
disp('I have computed the solution to the ODE by hand')
disp('However, I will use the built-in solver to improve my fluency')
x0=[0.1,0];
system_ode=@(t,x)[x1_dot(x(1),x(2),omega);x2_dot(x(1),x(2),omega)];
tspan=[0,5];
[t_ode, x_ode]=ode45(system_ode,tspan,x0);
figure;
plot(x_ode(:,1), x_ode(:,2), 'r-')
xlabel('x_1')
ylabel('x_2')
title('Phase Plane Plot using ode45')
grid on
axis equal
%Radius vs Time
syms r(t);
ode = diff(r,t) ==r-r^3
cond=r(0)==0.1
r_exact_sym=dsolve(ode,cond);
t_vec=linspace(0,5,500);
r_analytical=double(subs(r_exact_sym,t,t_vec));
r_numerical=sqrt(x_ode(:,1).^2+x_ode(:,2).^2);
figure;
plot(t_vec,r_analytical,'g-')
hold on;
plot(t_ode,r_numerical,'b--')
xlabel('Time t')
ylabel('Radius r')
title('Radius r(t) Convergence to Limit Cycle')
legend('Analytical Solution','Numerical (ode45)')
grid on;
hold off;
