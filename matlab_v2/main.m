% 20200418 by Dushan N. Wadduwage 
% added git
% main.m

clear
close all
% clc

%% read country data
[country_names, country_data] = read_country_data('world_data.js');
now_countryName              = 'Sri Lanka';%    'US';%          'Italy';%    'India'; %
population                   =  20E6;%          328.2E6;%       60.36E6;%     1.353E9;%

now_country_data             = country_data{find(contains(country_names,now_countryName)==1)}

start_date      = '2020-3-1';% '2020-1-22';%
start_date_ind  = find(contains({now_country_data.date},start_date)==1);
dates           = {now_country_data(start_date_ind(1):end).date};
gt.confirmed    = [now_country_data(start_date_ind(1):end).confirmed];
gt.deaths       = [now_country_data(start_date_ind(1):end).deaths]; 
gt.recovered    = [now_country_data(start_date_ind(1):end).recovered];
gt.active       = gt.confirmed - gt.deaths - gt.recovered;
gt.tvec         = 1:length(gt.active);

%%
S_0     = population;
T.E0    = 3;
T.E1    = 2;
T.I0    = 6;
T.I1    = 6;
T.I2    = 4;
T.I3    = 10;

fr.E0   = 1;
fr.E1   = 1;
fr.I0   = 0.3;
fr.I1   = 1 - fr.I0;
fr.R_I0 = 1;                   % 100% of I0 cases recover
fr.R_I1 = 0.8;                 % 80% of I1 cases recover
fr.I2   = 1 - fr.R_I1;         % rest goes to I2    
fr.R_I2 = 0.15/fr.I2;          % 15% of sever cases recover
fr.I3   = 1 - fr.R_I2;
fr.D    = 0.02/(fr.I3*fr.I2);  % 2% fatalify rate    
fr.R_I3 = 1 - fr.D;   

%%
rr.E0   = 0;
rr.E1   = 0;
rr.I0   = 0.001;
rr.I1   = 0.2;
rr.I2   = 1;
rr.I3   = 1;

t_evolve = gt.tvec(end) + 30;               % [days]
Beta    = [0.85*ones(1,13), 0.14*ones(1,21), 0.2*ones(1,7), 0.3*ones(1,t_evolve-41)] /S_0;
dt      = 1/24;

SL_population = POPULATION(S_0,T,fr,rr,Beta,dt)
SL_population.E0.N = 5;
evolutions = SL_population.evolve(t_evolve);

figure(1)
plot_results(evolutions, gt, t_evolve)

% error
%% optimizer

rng(2); %set random seed

% init 
model_population        = POPULATION(S_0,T,fr,rr,Beta,dt)
model_population.E0.N   = 5;

% optimization functon
options     = optimoptions('fmincon','MaxIterations',1000,'MaxFunctionEvaluations',2e6);
x0          = rand(1,length(gt.confirmed)+4);
x_lb        = [ zeros(1,length(gt.confirmed)) 0  0 1 1];
x_ub        = [2*ones(1,length(gt.confirmed)) .1 1 1 1];                       % x(end)=c and c E [0,1]

% wl regularizer
optF        = @(x) SEIR_objective_wlReg(x,model_population,gt)
tic
[x,fval]    = fmincon(optF,x0,[],[],[],[],x_lb,x_ub,[],options);
toc
[fval evolutions]  = optF(x);

%%
figure(2)
t_evolve = length(evolutions.S)
plot_results(evolutions, gt, t_evolve)

figure(3)
plot(x)



function plot_results(evolutions, gt, t_evolve)

FS = 14;

tvec = [1:t_evolve]';
pred_diagnosed = (evolutions.I0(2,:) + evolutions.I1(2,:) + evolutions.I2(2,:) + ...
                  evolutions.I3(2,:) + evolutions.R(2,:) + evolutions.D(2,:) )';

hold off
subplot(1,2,1);
plot(tvec, pred_diagnosed, 'k-x', gt.tvec, gt.confirmed, 'b-x')
set(gca,'fontsize',FS);
xlabel('Time [days]');
ylabel('No. of diagnosed cases');

legend('Predicted','Observed','location','NW');

subplot(1,2,2);
plot(tvec, sum(evolutions.D,1)', 'k-x', gt.tvec, gt.deaths, 'b-x')
set(gca,'fontsize',FS);
xlabel('Time [days]');
ylabel('No. of fatalities');

legend('Predicted','Observed','location','NW');


set(gcf, 'Units', 'Inches', 'Position', [1,1,14,6])

end





