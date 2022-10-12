clc,clear,close all

% Constants
xD = 7;      % Distance between turbine to rotor diameter fraction
% kw = 0.75;   % Wake expansion parameter (onshore)
kw = 0.04;   % Wake expansion parameter (offshore)

% parameter
a1 = linspace(0,1/3,1E5);  % Axial induction factor
a2 = a1;

% Power coeffcient
Cp_1 = 4.*a1.*(1-a1).^2;
Cp_2 = 4.*a2.*(1-a2).^2;

% Effciency function
eff = 1/2.*(Cp_1 + Cp_2.*(1 - (2*a1)./(1 + 2*kw*xD)^2).^3);

% Maxima of P --> a value
[eff_max,inx] = max(eff);

fprintf(['The maximum effciency of the system is %.3f,' ...
    ' having the axial induction factor of %.4f\n'],eff_max,a1(inx))

% Compare to operating at the Betz limit
[a_betz, eff_betz] = Betz_comparison(kw,xD);


% Effciency plot
figure()
hold on
plot(a1,eff,LineWidth=2)
max_p = plot(a1(inx),eff(inx),'.','MarkerSize',30);
plot(a_betz,eff_betz,'.','MarkerSize',30)

set(0,'defaultTextInterpreter','latex');
title('Effciency of two turbine farm for $a_1 = a_2$')
xlabel('Axial induction factor $a$')
ylabel('Farm effciency $\eta$')
legend('Effciency of wind farm','Maximum effciency','Betz limit','Interpreter','Latex','Location','northwest')
xlim([0 max(a1)+(0.35-1/3)])
ylim([0 max(eff)+0.1])
grid


function [a_betz,eff] = Betz_comparison(kw,xD)
% Determines the effciency of the system when operating at the Betz limit

% Axial induction factor at the Betz limit
a_betz = 1/3;   

% Power coeffcient
Cp_1 = 4.*a_betz.*(1-a_betz).^2;
Cp_2 = 4.*a_betz.*(1-a_betz).^2;

% Effciency function
eff = 1/2.*(Cp_1 + Cp_2.*(1 - (2*a_betz)./(1 + 2*kw*xD)^2).^3);

fprintf('The effciency of the system is %.3f, when Operating at the Betz limit\n',eff)

end

