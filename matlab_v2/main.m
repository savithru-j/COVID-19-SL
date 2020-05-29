% 20200418 by Dushan N. Wadduwage 
% added git
% main.m

clear all
close all
clc

%% read country data
[country_names, country_data] = read_country_data('../world_data.js');
now_countryName              = 'Sri Lanka';%    'US';%          'Italy';%    'India'; %
population                   =  21323733;%          328.2E6;%       60.36E6;%     1.353E9;%

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
T.E0    = 3;
T.E1    = 2;
T.I0    = 6;
T.I1    = 6;
T.I2    = 4;
T.I3    = 10;

fr.E0   = 1;
fr.E1   = 1;
fr.I0   = 0.30;
fr.I1   = 1 - fr.I0;
fr.R_I0 = 1;                   % 100% of I0 cases recover
fr.R_I1 = 0.80;                % 80% of I1 cases recover
fr.I2   = 1 - fr.R_I1;         % rest goes to I2    
fr.R_I2 = 0.75;                % 75% of I2 cases recover (1 - 0.8)*0.75 = 0.15
fr.I3   = 1 - fr.R_I2;         % rest goes to I3
fr.D    = 0.02/(fr.I2*fr.I3);  % 2% case fatality rate    
fr.R_I3 = 1 - fr.D;
% fr.R_I3 = 0.60;                % 60% of I3 cases recover (1 - 0.8)*(1 - 0.75)*0.6 = 0.03
% fr.D    = 1 - fr.R_I3;         % rest is fatal: (1 - 0.8)*(1 - 0.75)*(1 - 0.6) = 0.02

%%
S_0     = population;
t_evolve= gt.tvec(end) + 14;               % [days]
Beta    = zeros(t_evolve,1);
Beta(:) = 0.3/S_0;                         % legacy value: [0.85*ones(1,13), 0.14*ones(1,21), 0.2*ones(1,7), 0.3*ones(1,t_evolve-41)] /S_0; 
c       = zeros(t_evolve,4);
c(:,1)  = 0;                               % c_I0(t) 
c(:,2)  = .1;                              % c_I1(t) 
c(:,3)  = 1;                               % c_I2(t) 
c(:,4)  = 1;                               % c_I3(t) 
dt      = 1/24;

SL_population = POPULATION(S_0,T,fr,Beta,c,dt)
SL_population.E0.N = 5;
SL_population.R.N(2) = 1;
SL_population.S.N = S_0 - SL_population.E0.N - SL_population.R.N(2);
evolutions = SL_population.evolve(t_evolve);
% 
figure(1)
plot_results(evolutions, gt, t_evolve,dates)

%% optimizer
rng(2);                                                     %set random seed
model_population        = POPULATION(S_0,T,fr,Beta,c,dt)    % init 
model_population.E0.N   = 5;
model_population.R.N(2) = 1;
model_population.S.N = model_population.S_0 - model_population.E0.N - model_population.R.N(2);

options     = optimoptions('fmincon',...
                            'MaxIterations',10000,...
                            'MaxFunctionEvaluations',1e6,...
                            'FiniteDifferenceStepSize',0.001,...
                            'UseParallel',false);
x0          = rand(length(gt.confirmed),6);
x_lb        = zeros(size(x0));
x_ub        = zeros(size(x0));
x_lb(:,1)   = 0;            % transmission rate
x_ub(:,1)   = 2; 
x_lb(:,2)   = 0;            % c_I0
x_ub(:,2)   = 1;
x_lb(:,3)   = 0;            % c_I1
x_ub(:,3)   = 1;
x_lb(:,4)   = 1;            % c_I2
x_ub(:,4)   = 1;
x_lb(:,5)   = 1;            % c_I3
x_ub(:,5)   = 1;

%                                                                    case fatality ratio
%                                                          recovery fraction from I2
%                                                  recovery fraction from I1
%                                            asymptomatic fraction   
%              <------time durations-------->                 

% index        1    2    3    4    5    6    7     8       9        10 
% typical      3    2    6    6    4    10   0.3   0.8     0.75     0.60 
% paramter     E0.T E1.T I0.T I1.T I2.T I3.T fr.I0 fr.R_I1 fr.R_I2  CFR
x_lb(1:10,6)= [3    2    6    6    4    10   0.3   0.1     0.1      0.0];
x_ub(1:10,6)= [3    2    30   30   4    10   0.3   0.9     0.9      0.02];

% x_lb(1:10,6)= [ 1    1    1    1    1    1   0.3   0.1     0.1      0.0];
% x_ub(1:10,6)= [10   10   10   10   10   10   0.3   0.9     0.9      0.02];


err_type    = 'L2_type1';   % {'log_type1','log_type2','L2_type1','L2_type2'}

% wl regularizer


optF        = @(x) SEIR_objective_wlReg(x,model_population,gt,err_type)  % optimization functon
parfor i=1:40
    i
    tic
    x0          = x_lb + (x_ub-x_lb).*rand(length(gt.confirmed),6);
    [xx{i},fval]    = fmincon(optF,x0(:),[],[],[],[],x_lb(:),x_ub(:),[],options);
    toc
end
x = xx{29};


[fval, evolutions]  = optF(x);

%% Figures
figure(2)
t_evolve = length(evolutions.S)-1
plot_results(evolutions, gt, t_evolve,dates)
%saveas(gca,sprintf('./_Figs/fit.jpeg'));

figure(3)
FS = 14;
x = reshape(x,t_evolve,6);
x(1:10,6)
plot(x(:,1:5),'linewidth',2)
set(gca,'fontsize',FS);
xlabel('Time [days]');
ylabel('Rate [AU]');
xticks(1:t_evolve)
xticklabels(dates)
xtickangle(90)
legend('Beta_1','c_{I0}','c_{I1}','c_{I2}','c_{I3}');
set(gcf, 'Units', 'Inches', 'Position', [1,1,14,6])
%saveas(gca,sprintf('./_Figs/paramteres.jpeg'));



