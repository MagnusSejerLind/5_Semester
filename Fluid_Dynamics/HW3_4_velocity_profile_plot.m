clc,clear,close all


% Position between plates(y/a)
ya = linspace(0,1);

% Dimensionless velocity profiles
up_top = 2.*ya - ya.^2;
up_bot =  ya.^2;


f = figure();
f.Position = [295 250 900 500];
set(0,'defaultTextInterpreter','latex'); 
hold on

plot(up_top,ya,'b',LineWidth=2)
plot(up_bot,ya,'r',LineWidth=2)

title('Dimensionless velocity of laminar flow between plates with the upper plate moving at constant velocity')
legend('Upper plate','Lower plate','Interpreter','latex',Location='northwest');
xlabel('$\frac{u}{U}$')
ylabel('$\frac{y}{a}$',Rotation=360)
grid





