
%% load data %% 

clearvars

ix = 1;
tf = dir('*time_*');
t = load(tf(ix).name);

M2f = dir('*M2_*'); 
M2f(ix).date
M2 = load(M2f(ix).name);

M4f = dir('*M4_*'); 
M4f(ix).date
M4 = load(M4f(ix).name);

Z13f = dir('*Z13_*');Z13f(ix).date
Z13 = load(Z13f(ix).name);

Z22f = dir('*Z22_*'); 
ix = 1;
Z22f(ix).date
Z22 = load(Z22f(ix).name);

%% simulation vars %%

abar = 6; 
D0=1;
P = .4; 
zeta = 1/abar/P;
tauR = abar * 0.5 / P;
Dinf = D0/(1+zeta);

Dt = M2/2./t-Dinf; 
kurt = M4./M2.^2-3;

%% fit specifications %%

theta = 1; 

%time-point bounds (seconds)
t_lb = 10;
t_ub = 100;

%index range for fit
t1 = find(t <= t_lb); 
n1 = t1(end);           %index for lower bound of data for fit
t2 = find(t >= t_ub); 
n2 = t2(1);             %index for upper bound of data for fit

%plot 1 fit
xdata1 = t(n1:n2)/tauR; 
ydata1 = Dt(n1:n2); 
fun1 = @(f,xdata1)f(1) + f(2)*xdata1.^(-theta); 
f0 = [1,1]; 
f = lsqcurvefit(fun1, f0, xdata1, ydata1); 

%plot 2 fit
xdata2 = t(n1:n2)/tauR; 
ydata2 = kurt(n1:n2); 
fun2 = @(g,xdata2)g(1) + g(2)*xdata2.^(-theta); 
g0 = [1,1]; 
g = lsqcurvefit(fun2, g0, xdata2, ydata2); 

%plot 3 fit
xdata3 = t(n1:n2)/tauR; 
ydata3 = Z13(n1:n2); 
fun3 = @(h,xdata3)h(1) + h(2)*xdata3.^(-theta); 
h0 =[1,1]; 
h = lsqcurvefit(fun3, h0, xdata3, ydata3); 

%plot 4 fit
xdata4 = t(n1:n2)/tauR; 
ydata4 = Z22(n1:n2); 
fun4 = @(j,xdata4)j(1) + j(2)*xdata4.^(-theta); 
j0 =[1,1]; 
j = lsqcurvefit(fun4, j0, xdata4, ydata4); 

%% figures %%

%figure 1
subplot(2,2,1)
loglog(t/tauR, abs(Dt),'.') 
hold on
fplot(@(t)(f(1) + f(2)*t.^(-theta)),[t(n1),t(n2)])
title('D(t)')
ylabel('D(t)-D_\infty')
xlabel('t/\tau_R')

%figure 2
subplot(2,2,2)
loglog(t(kurt>0)/tauR, kurt(kurt>0),'.')
hold on
fplot(@(t)(g(1) + g(2)*t.^(-theta)),[t(n1),t(n2)])
title('K(t)')
ylabel('K(t)')
xlabel('t/\tau_R')

%figure 3
subplot(2,2,3)
plot(t(1:length(Z13))/tauR, Z13,'.')
hold on 
fplot(@(t)(h(1) + h(2)*t.^(-theta)),[t(n1),t(n2)])
title('Z^{(1,3)}(t)')
ylabel('Z^{(1,3)}(t)')
xlabel('t')

%figure 4
subplot(2,2,4)
plot(t(1:length(Z22))/tauR, Z22,'.')
hold on 
fplot(@(t)(j(1) + j(2)*t.^(-theta)),[t(n1),t(n2)])
title('Z^{(2,2)}(t)')
ylabel('Z^{(2,2)}(t)')
xlabel('t')
