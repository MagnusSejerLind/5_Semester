%% Course project - Renewable Energy Networks
clc,clear,close all

plot_option = 1;    % Set to 1 to plot


% Import data
[imp_data] = import_data;

% Setup data
[data] = setup_data(imp_data);
clear imp_data

% Nomalize data
[data] = Normalize_data(data);

% Determine mismatch
[data] = Mismatch(data,plot_option);

% Pricipal component analysis
[eigen] = PCA(data,plot_option);

% Plot eigenvectors
PCAplot(eigen,plot_option)




% Import data
function [imp_data] = import_data
imp_data.offshore = importdata('Data\offshore_wind_1979-2017.csv');
imp_data.onshore = importdata('Data\onshore_wind_1979-2017.csv');
imp_data.elec_demand = importdata('Data\electricity_demand.csv');
imp_data.heat_demand = importdata('Data\heat_demand.csv');
imp_data.pv_CF = importdata('Data\pv_optimal.csv');
end




function [data] = setup_data(imp_data)
 
% Only use onshore, electricity & PV to be used

%% Offshore CF

% Time, Hourly UTC time, 01/01-1979,00:00 to 32/12-2017,23:00
data.offshore.time = imp_data.offshore.textdata(2:end,1);   

% All CF
data.offshore.CF_all = imp_data.offshore.data;              

% Country specific CF 
fn_in_country_offshore = imp_data.offshore.textdata(1,2:end);
for i = 1:length(fn_in_country_offshore)
    data.offshore.country.(fn_in_country_offshore{i}).CF = imp_data.offshore.data(1:end,i);
end


%% Onshore CF

% Time, Hourly UTC time, 01/01-1979,00:00 to 32/12-2017,23:00
data.onshore.time = imp_data.onshore.textdata(2:end,1);

% All CF
data.onshore.CF_all = imp_data.onshore.data;           

% Country specific CF 
fn_in_country_onshore = imp_data.onshore.textdata(1,2:end);
for i = 1:length(fn_in_country_onshore)
    data.onshore.country.(fn_in_country_onshore{i}).CF = imp_data.onshore.data(1:end,i);
end


%% PV CF

% Time, Hourly UTC time, 01/01-1979,00:00 to 32/12-2017,23:00
data.pv.time = imp_data.pv_CF.textdata(2:end,1);

% All CF
data.pv.CF_all = imp_data.pv_CF.data;

% Country specific CF 
fn_in_country_pv = imp_data.pv_CF.textdata(1,2:end);
for i = 1:length(fn_in_country_pv)
    data.pv.country.(fn_in_country_pv{i}).CF = imp_data.pv_CF.data(1:end,i);
end


%% Heating demand

% Time, Hourly UTC time, 01/01-2015,00:00 to 32/12-2015,23:00
data.heat.time = imp_data.heat_demand.textdata(2:end,1);

% All demand
data.heat.demand_all = imp_data.heat_demand.data;

% Country specific demand
fn_in_country_heat = imp_data.heat_demand.textdata(1,2:end);
for i = 1:length(fn_in_country_heat)
    data.heat.country.(fn_in_country_heat{i}).demand = imp_data.heat_demand.data(1:end,i);
end

%% Electricity demand


% Time, Hourly UTC time, 01/01-2015,00:00 to 32/12-2015,23:00
data.elec.time = imp_data.elec_demand.textdata(2:end,1);

% All demand
data.elec.demand_all = imp_data.elec_demand.data;

% Country specific demand
fn_in_country_elec = imp_data.elec_demand.textdata(1,2:end);
for i = 1:length(fn_in_country_elec)
    data.elec.country.(fn_in_country_elec{i}).demand = imp_data.elec_demand.data(1:end,i);
end


%% Modify data form non-represented countries

% Some value are unkown, these are set to zero  

% All values for heat and pv are known, the specific missing countries for the
% other sectors are determined
unknown.offshore = setxor(fn_in_country_heat, fn_in_country_offshore);
unknown.onshore = setxor(fn_in_country_heat, fn_in_country_onshore);
unknown.elec = setxor(fn_in_country_heat, fn_in_country_elec);


% Offshore
for i = 1:length(unknown.offshore)
    data.offshore.country.(unknown.offshore{i}).CF = zeros(length(data.offshore.time),1);
end

% Onshore
for i = 1:length(unknown.onshore)
    data.onshore.country.(unknown.onshore{i}).CF = zeros(length(data.onshore.time),1);
end

% Electricity
data.elec.country.(unknown.elec{:}).demand = zeros(length(data.elec.time),1);


% Sort countries alphabetical
data.elec.country = orderfields(data.elec.country);
data.onshore.country = orderfields(data.onshore.country);
data.pv.country = orderfields(data.pv.country);




%% Extract single year - 2015

fn_in_country = fieldnames(data.onshore.country);

for i = 1:length(fn_in_country)
    k = 0;
    for j = 315577:315577 + (8760 - 1)
    % 315577 index for 1/1/2015
    k = k + 1;
        
        % Onshore
        data.onshore.country.(fn_in_country{i}).CF_2015(k,1) = data.onshore.country.(fn_in_country{i}).CF(j);
        
        % PV
        data.pv.country.(fn_in_country{i}).CF_2015(k,1) = data.pv.country.(fn_in_country{i}).CF(j);
    end
end

%% Generation 2015 
% CF * Demand

for i = 1:length(fn_in_country)
    % Wind generation
    data.onshore.country.(fn_in_country{i}).gen_2015 = data.onshore.country.(fn_in_country{i}).CF_2015 .* data.elec.country.(fn_in_country{i}).demand;

    % PV generation
    data.pv.country.(fn_in_country{i}).gen_2015 = data.pv.country.(fn_in_country{i}).CF_2015 .* data.elec.country.(fn_in_country{i}).demand;
end



end








function [data] = Normalize_data(data)
% Normalizes the wind and PV 2015 generation to the load

fn_country = fieldnames(data.elec.country);

for i = 1:length(fn_country)
% x = (x*mean(y))/mean(x)

    % Wind
    data.onshore.country.(fn_country{i}).norm_gen_2015 = ( data.onshore.country.(fn_country{i}).gen_2015 .* mean(data.elec.country.(fn_country{i}).demand)  ) ./ mean(data.onshore.country.(fn_country{i}).gen_2015);

    % PV
    data.pv.country.(fn_country{i}).norm_gen_2015 = ( data.pv.country.(fn_country{i}).gen_2015 .* mean(data.elec.country.(fn_country{i}).demand)  ) ./ mean(data.pv.country.(fn_country{i}).gen_2015);


% x = (x - min(x)) / (max(x) - min(x))
%     data.elec.country.(fn_coun{i}).norm_demand = ( data.elec.country.(fn_coun{i}).demand - min(data.elec.country.(fn_coun{i}).demand) ) / ( max(data.elec.country.(fn_coun{i}).demand) - min(data.elec.country.(fn_coun{i}).demand));

end

end

function [data] = Mismatch(data,plot_option)
% Determines the mismatch of 2015
clc

alpha = (0:0.2:1);


fn_alpha = ["alpha_00" 'alpha_02' 'alpha_04' 'alpha_06' 'alpha_08' 'alpha_1'];


fn_country = fieldnames(data.onshore.country);
for i = 1:length(fn_country)
    for j = 1:length(fn_alpha)

        % (alpha*G_w) + ((1-alpha)*G_s) - L
        data.mismatch_2015.(fn_country{i}).(fn_alpha{j}) = ( alpha(j)*data.onshore.country.(fn_country{i}).norm_gen_2015 ) + (1 - alpha(j)) .* data.pv.country.(fn_country{i}).norm_gen_2015 - data.elec.country.(fn_country{i}).demand;

        % Set NAN to 0
        data.mismatch_2015.(fn_country{i}).(fn_alpha{j})(isnan(data.mismatch_2015.(fn_country{i}).(fn_alpha{j}))) = 0;
    end
end




%% Plot mismatch for DEU & DNK

if plot_option == 1

fn_title = ["$\alpha = 0.0 $" '$\alpha = 0.2$' '$\alpha = 0.4 $' '$\alpha = 0.6 $' '$\alpha = 0.8 $' '$\alpha = 1.0 $'];

% DEU
figure()
tl_DEU = tiledlayout(2,3);
title(tl_DEU,'Mismatch for Germany','Interpreter','latex')

for i = 1:length(fn_alpha)
    nexttile
    plot(data.mismatch_2015.DEU.(fn_alpha{i}))
    
    title(fn_title{i},'Interpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    xlabel('Time [h]')
    ylabel('$\Delta$',Rotation=360)
    xlim([0 8760]);
    ylim([-100000 max(data.mismatch_2015.DEU.alpha_00)+5E4])
    grid
end
print('Mismatch_DEU', '-depsc');  


% DNK
figure()
tl_DNK = tiledlayout(2,3);
title(tl_DNK,'Mismatch for Denmark','Interpreter','latex')
for i = 1:length(fn_alpha)

    nexttile
    plot(data.mismatch_2015.DNK.(fn_alpha{i}))
  
    title(fn_title{i},'Interpreter','latex');
    xlabel('Time [h]')
    ylabel('$\Delta$',Rotation=360)
    xlim([0 8760]);
    ylim([-7E3 2.5E4])
    grid
end
print('Mismatch_DNK', '-depsc');  

end

%% Aggregated European mismatch



for i = 1:length(fn_alpha)
    data.EU_mismatch_2015.(fn_alpha{i}) = 0;
    for j = 1:length(fn_country)
        country_mismatch_2015 = data.mismatch_2015.(fn_country{j}).(fn_alpha{i});

        data.EU_mismatch_2015.(fn_alpha{i}) = data.EU_mismatch_2015.(fn_alpha{i}) + country_mismatch_2015;

    end
end


%% Plot aggregated European mismatch

if plot_option == 1
figure()
tl_EU = tiledlayout(2,3);
title(tl_EU,'Aggregated mismatch for Europe','Interpreter','latex')

for i = 1:length(fn_alpha)
    nexttile
    plot(data.EU_mismatch_2015.(fn_alpha{i}))
    
    title(fn_title{i},'Interpreter','latex');
    set(0,'defaultTextInterpreter','latex');
    xlabel('Time [h]')
    ylabel('$\Delta$',Rotation=360)
    xlim([0 8760]);
    ylim([min(data.EU_mismatch_2015.alpha_00)-5E4 max(data.EU_mismatch_2015.alpha_00)+5E4])
    grid
end
print('Mismatch_EU', '-depsc');  

end



%% Mismatch variance

variance_DNK_mismatch_2015 = zeros(length(fn_alpha),1);
variance_DEU_mismatch_2015 = zeros(length(fn_alpha),1);
variance_EU_mismatch_2015 = zeros(length(fn_alpha),1);


% Determine the variance of Denmark and Germany
for i = 1:length(fn_alpha)
    variance_DNK_mismatch_2015(i) = var(data.mismatch_2015.DNK.(fn_alpha{i}));
    variance_DEU_mismatch_2015(i) = var(data.mismatch_2015.DEU.(fn_alpha{i}));
    variance_EU_mismatch_2015(i) = var(data.EU_mismatch_2015.(fn_alpha{i}));

end

% Variance with respect to mean load^2
var_load_DNK = variance_DNK_mismatch_2015 / mean(data.elec.country.DNK.demand)^2;
var_load_DEU = variance_DEU_mismatch_2015 / mean(data.elec.country.DEU.demand)^2;
var_load_EU = variance_EU_mismatch_2015 / mean(mean(data.elec.demand_all)).^2;


if plot_option == 1

figure()
plot(alpha,var_load_DNK,'.-',MarkerSize=20,LineWidth=2)

title('Mismatch Variance Denmark')
xlabel('$\alpha$')
ylabel('Var$(\vec{\Delta}) / \langle L_{\mbox{DNK}}\rangle ^2$')
ylim([0 max(var_load_DNK)+0.5])
grid
print('Mismatch_variance_DNK', '-depsc');  



figure()
plot(alpha,var_load_DEU,'.-',MarkerSize=20,LineWidth=2)

title('Mismatch Variance Germany')
xlabel('$\alpha$')
ylabel('Var$(\vec{\Delta}) / \langle L_{\mbox{DNK}}\rangle ^2$')
ylim([0 max(var_load_DEU)+0.5])
grid
print('Mismatch_variance_DEU', '-depsc');  


figure()
plot(alpha,var_load_DEU,'.-',MarkerSize=20,LineWidth=2)

title('Mismatch Variance aggregated Europe')
xlabel('$\alpha$')
ylabel('Var$(\vec{\Delta}) / \langle L_{\mbox{DNK}}\rangle ^2$')
grid
print('Mismatch_variance_EU', '-depsc');  



end


end






function [eigen] = PCA(data,plot_option)
% Principal component analysis


fn_country = fieldnames(data.mismatch_2015);
fn_alpha = ["alpha_00" 'alpha_02' 'alpha_04' 'alpha_06' 'alpha_08' 'alpha_1'];
alpha = (0:0.2:1);


% Mismatch
for i = 1:length(fn_country)
    % 32 countries x 8760 mismatch h data
    mm_a0(:,i) = data.mismatch_2015.(fn_country{i}).alpha_00;
    mm_a02(:,i) = data.mismatch_2015.(fn_country{i}).alpha_02;
    mm_a04(:,i) = data.mismatch_2015.(fn_country{i}).alpha_04;
    mm_a06(:,i) = data.mismatch_2015.(fn_country{i}).alpha_06;
    mm_a08(:,i) = data.mismatch_2015.(fn_country{i}).alpha_08;
    mm_a1(:,i) = data.mismatch_2015.(fn_country{i}).alpha_1;

end

% Determine covarince matrices
mismatch_cov.(fn_alpha{1}) = cov(mm_a0);
mismatch_cov.(fn_alpha{2})= cov(mm_a02);
mismatch_cov.(fn_alpha{3}) = cov(mm_a04);
mismatch_cov.(fn_alpha{4}) = cov(mm_a06);
mismatch_cov.(fn_alpha{5}) = cov(mm_a08);
mismatch_cov.(fn_alpha{6}) = cov(mm_a1);

% Determine eigenvalues and eigenvectors
for i = 1:6
    eigen.vec.(fn_alpha{i}) = 0;
    eigen.val_mat.(fn_alpha{i}) = 0;
    [eigen.vec.(fn_alpha{i}),eigen.val_mat.(fn_alpha{i})] = eig(mismatch_cov.(fn_alpha{i}));

    % Extract eigenvalues
    eigen.val.(fn_alpha{i}) = diag(eigen.val_mat.(fn_alpha{i}));

    % Sort in decending order
    eigen.val.(fn_alpha{i}) = flip(eigen.val.(fn_alpha{i}));
    eigen.vec.(fn_alpha{i}) = flip(eigen.vec.(fn_alpha{i}),2);

    % Cumulative sum of eigenvalues
    eigen.cum_val.(fn_alpha{i}) = cumsum(eigen.val.(fn_alpha{i}));

    % Set sum to 1
    eigen.cum_val_1.(fn_alpha{i}) = eigen.cum_val.(fn_alpha{i}) / max(eigen.cum_val.(fn_alpha{i}));
    eigen.val_1.(fn_alpha{i}) = eigen.val.(fn_alpha{i}) / max(eigen.val.(fn_alpha{1}));

end



eigen_val_mat = zeros(length(fn_alpha),length(fn_alpha));
% eigen_val_mat_scale = zeros(length(fn_alpha),length(fn_alpha));
for i = 1:length(fn_alpha)
    % Create matrix with the six largest eigenvalues for each alpha increment
    eigen_val_mat(:,i) = eigen.val.(fn_alpha{i})(1:6);

    % Cumulative eigenvalues
    eigen_val_mat_cum_unscale(:,i) = eigen.cum_val.(fn_alpha{i})(1:6);
    eigen_val_mat_cum_tot(:,i) = eigen.cum_val.(fn_alpha{i});

    % Cumulative eigenvalues scaled
    eigen_val_mat_cum(:,i) = eigen_val_mat_cum_unscale(:,i) ./ max(eigen_val_mat_cum_tot(:,i));

end


if plot_option == 1


    % Determine the change between principal components
    mat(:,1) = eigen_val_mat_cum(1,:);
    mat(:,2) = eigen_val_mat_cum(2,:) -  eigen_val_mat_cum(1,:);
    mat(:,3) = eigen_val_mat_cum(3,:) - eigen_val_mat_cum(2,:);
    mat(:,4) = eigen_val_mat_cum(4,:) - eigen_val_mat_cum(3,:);
    mat(:,5) = eigen_val_mat_cum(5,:) - eigen_val_mat_cum(4,:);
    mat(:,6) = eigen_val_mat_cum(6,:) - eigen_val_mat_cum(5,:);
    mat(:,7) = 1 - eigen_val_mat_cum(6,:);

    % impact of eigenvalues on system
    figure()
    area(alpha,mat)

    title('PCA eigenvalues')
    xlabel('$\alpha$')
    ylabel('$\lambda_k$')
    legend('$k=1$','$k=2$','$k=3$','$k=4$','$k=5$','$k=6$','$7 \leq k \leq 32$','Interpreter','latex',Location='best')
    print('PCA_eigenvalue', '-depsc');


end

end




function PCAplot(eigen,plot_option)  

if plot_option == 1

% country names
cou_nam = ({'Belgium', 'Bulgaria', 'Czech Republic', 'Denmark', 'Germany','Estonia','Ireland',...
    'Greece','Spain','France','Croatia','Italy', 'Cyprus','Latvia', 'Lithuania', 'Luxembourg',...
     'Hungary', 'Malta', 'Netherlands','Austria', 'Poland', 'Portugal','Romania','Slovenia', 'Slovakia',...
        'Finland', 'Sweden', 'Bosnia and Herzegovina', 'Malta'});
numRegions = length(cou_nam);

PC_nam = ["$\lambda_1$" "$\lambda_2$" "$\lambda_3$" "$\lambda_4$" "$\lambda_5$" "$\lambda_6$"];
set(0,'defaultTextInterpreter','latex');
%% Alpha = 0.0

figure()
tl_alpha_00 = tiledlayout(1,3);
% title(tl_alpha_00,'Top Principal components for $\alpha=0.0$')

for j = 1:3
% Plots the first three principal components    
nexttile

% Set color value for range [-1 1]
color_val_alpha_00 = round(interp1([-1 1],[1 30],eigen.vec.alpha_00(:,j)));

worldmap("Europe");
geoshow('landareas.shp');
bordersm('countries','facecolor','white');
cmap = turbo(30); % Number of color increments
colormap(cmap)


for i = 1:numRegions
        bordersm(cou_nam{i},'facecolor',cmap(color_val_alpha_00(i),:)); 
end

title(PC_nam(j));

end

colorbar
set(gcf, 'Position', get(0, 'Screensize'));
print('Principal_components_alpha_00', '-depsc');  


%% Alpha = 1.0

figure()
tl_alpha_1 = tiledlayout(2,3);
% title(tl_alpha_1,'Top Principal components for $\alpha=1.0$')

for j = 1:6
% Plots the first three principal components
nexttile

% Set color value for range [-1 1]
color_val_alpha_1 = round(interp1([-1 1],[1 30],eigen.vec.alpha_1(:,j)));

worldmap("Europe");
geoshow('landareas.shp');
bordersm('countries','facecolor','white');
cmap = turbo(30); 

for i = 1:numRegions
        bordersm(cou_nam{i},'facecolor',cmap(color_val_alpha_1(i),:)); 
end

colormap(cmap)
title(PC_nam(j));

end

colorbar
set(gcf, 'Position', get(0, 'Screensize'));
print('Principal_components_alpha_1', '-depsc');  

end
end


