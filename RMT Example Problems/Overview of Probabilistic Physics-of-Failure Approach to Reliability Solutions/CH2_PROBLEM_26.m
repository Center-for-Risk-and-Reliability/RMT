% CHAPTER 2 PROBLEM 26 - SCC
% Script by Dr. Mehdi Amiri
%========Stress Corrosion Cracking Modeling===========
% This code will simulate the results that are provided in section
% 4.5.2 "Stress Corrosion Cracking (SCC) Modeling" of the final report.
% The model is based on Shoji's approach. Model description is given in
% section 4.5.2.
%============================================

%Model parameters for Aluminum 7075-T6 for data from Campbell 1982
M = 27; n = 1.53; m = 0.4; r0 = 4e-5; i0 = 0.0005; t0 = 0.1; epsf = 1e-3;
beta = 5.81; lambda = 0.2; z = 3; F = 96500; rhu = 2.7e3; E = 71e3;
sigmy = 503;

Kdot = 0;
Coeff1 = (M*i0/(z*rhu*F*(1-m)))*(t0/epsf)^m;
Coeff11 = (M*i0/(100*z*rhu*F*(1-m)))*(t0/epsf)^m;
Coeff2 = beta*sigmy*n/((r0^m)*E*(n-1));
Coeff22 = beta*sigmy*n/(r0*E*(n-1));
Ki = 7;                  %Starting K value for simulations
Kf = 30;                 %Final K value for simulations
K = [Ki:Kf];             %Range of K values for simulations

count = 1;
for i = Ki:Kf
dadt(count) = 1e-3*(Coeff11*(Coeff22*(log(lambda*((i/sigmy)^2)/r0))^(1/(n-1))))^(m/(1-m));
    count = count + 1;
end

%=========Plot modeling results==============
figure(1)
p1 = semilogy(K,dadt,'LineWidth',2);
axis([0 30 1.00E-11 1.00E-07]);
set(gca,'fontsize',14);
grid on
xlabel('\itK\rm [MPa\cdotm^{1/2}]')
ylabel('d\ita\rm/d\itt\rm [m/s]'); 
hold on
%====================================================

%experimental data for SCC of alloy 2014, wet twice a day with 3.5% NaCI, 25”C
expdata =[5.7777777	1.17E-09
6.5555553	2.56E-09
8.944445	1.17E-09
8.944445	1.90E-09
8.944445	2.72E-09
9           4.40E-09
10.388889	1.82E-09
10.388889	3.53E-09
12.222222	3.46E-09
12.111111	5.27E-09
12.944445	4.23E-09
13          5.06E-09
17.166666	5.06E-09
17.222221	1.13E-08
18.666666	5.38E-09
24.11111	1.13E-08];


%===========Plot experimental data=======================
p2 = semilogy(expdata(:,1),expdata(:,2), 'r*', 'MarkerSize',12);
 legend([p1 p2], 'SCC model prediction', 'Exp. Al-alloy 2014 in 3.5% NaCl at 25^oC');