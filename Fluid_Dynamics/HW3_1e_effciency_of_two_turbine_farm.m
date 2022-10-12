clc,clear,close all

% Constants
xD = 7;      % Distance between turbine to rotor diameter fraction
kw = 0.075;   % Wake expansion parameter (onshore)
% kw = 0.04;   % Wake expansion parameter (offshore)

% parameter
a1 = linspace(0,1/3,100);  % Axial induction factor
a2 = a1;

eff = zeros(length(a1));
eff_max = 0;

for i = 1:length(a1)
    for j = 1:length(a2)

        % Power coeffcient
        Cp_1 = 4.*a1(i).*(1-a1(i)).^2;
        Cp_2 = 4.*a2(j).*(1-a2(j)).^2;

        % Effciency function
        eff(i,j) = 1/2.*(Cp_1 + Cp_2.*(1 - (2*a1(i))./(1 + 2*kw*xD)^2).^3);

        if max(eff(:)) > eff_max
            eff_max = max(eff(:));
            a1_max = a1(i);
            a2_max = a2(j);
            a1_max_ind = i;
            a2_max_ind = j;
        end
    end
end



fprintf(['The maximum effciency of the system is %.3f,' ...
    ' having the axial induction factor of a1 = %.4f, a2 = %.4f\n'],eff_max,a1_max,a2_max)


% Compare to operating at the Betz limit
[a_betz, eff_betz] = Betz_comparison(kw,xD);


% Figure
f = figure();
sc = 300; % Axis scaling
f.Position = [200 250 1000 500];

s = surf(eff); % Effciency surface plot
hold on
plot3(a2_max_ind,a1_max_ind,eff_max,'k.','MarkerSize',50) % Maximum effciency
plot3(a_betz*sc,a_betz*sc,eff_betz,'m.','MarkerSize',50) % Betz limit comparison

set(0,'defaultTextInterpreter','latex');
legend('Effciency','Maximum effciency','Effciency at the Betz limit','Location','northeast')
title('Effciency of two turbine farm')
s.EdgeColor = 'none';
colorbar
ylabel('$a_1$')
xlabel('$a_2$')
zlabel('Effciency')
xticks([0 20 40 60 80 100])
xticklabels({'0','0.066','0.133','0.200','0.266','0.333'})
yticks([0 20 40 60 80 100])
yticklabels({'0','0.066','0.133','0.200','0.266','0.333'})




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

