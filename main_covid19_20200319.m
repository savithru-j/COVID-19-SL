% paramters from references (https://alhill.shinyapps.io/COVID19seir/)
clc; clear all; close all

IncubPeriod     = 5;                % 1/a
DurMildInf      = 6;                % 1/(p1+g1)
DurSevere       = 4;                % 1/(p2+g2)     [Time from symptoms to ICU admit=10] - [Duration of mild infections=6]
DurICU          = 10;               % 1/(u+g3)      [Time from hospital admit to death=14] - [Duration of severe infections=4]

% known probabilities from lit
prob_I1_E       = 1;                %               [All exposed develop to mild cases]    
prob_R_I1       = 0.81*prob_I1_E;   % g1/(p1+g1)    [81% of the cases are mild]
prob_I2_I1      = 1 - prob_R_I1;    % p1/(p1+g1)    
prob_R_I2       = 0.14/prob_I2_I1;  % g2/(p2+g2)    [14% of the cases are severe]
prob_I3_I2      = 1 - prob_R_I2;    % p2/(p2+g2)
% prob_I3_I2      = 0.06/prob_I2_I1;                [6% of the cases are critical], Note: This number doesnt addup. Look at lterature 
prob_D_I3       = 0.02/(prob_I3_I2*prob_I2_I1);  % u/(p2+g2)     [2% die after being crtical, assuming after being in the ICU]
prob_R_I3       = 1-prob_D_I3;      % u/(p2+g2)

CFR             = 0.02;             %               [Case fatality ratio]             

% rates
a   = (1/IncubPeriod)*prob_I1_E;    % [day^-1]
g1  = (1/DurMildInf) * prob_R_I1;
p1  = (1/DurMildInf) * prob_I2_I1;
g2  = (1/DurSevere)  * prob_R_I2;
p2  = (1/DurSevere)  * prob_I3_I2;
g3  = (1/DurICU)     * prob_R_I3;
u   = (1/DurICU)     * prob_D_I3;

% transmission rates. Beta values are always scaled by N
b1N = 1;                            % [day^-1]     bi rate at which infected individuals in class Ii contact susceptibles and infect them
b2N = 0;
b3N = 0;

N   = 1.3e9;                                         % population

b1  = b1N/N; 
b2  = b2N/N;
b3  = b3N/N;

% equations
% N   = S+I1+I2+I3+R+D;                              % population doesnt change 
dS  = @(S,E,I1,I2,I3,R,D) -b1*I1*S - b2*I2*S - b3*I3*S; % dS = dS/dt   
dE  = @(S,E,I1,I2,I3,R,D) b1*I1*S + b2*I2*S + b3*I3*S - a*E;           
dI1 = @(S,E,I1,I2,I3,R,D) a*E - g1*I1 - p1*I1;
dI2 = @(S,E,I1,I2,I3,R,D) p1*I1 - g2*I2 - p2*I2;
dI3 = @(S,E,I1,I2,I3,R,D) p2*I2 - g3*I3 - u*I3;
dR  = @(S,E,I1,I2,I3,R,D) g1*I1 + g2*I2 + g3*I3;
dD  = @(S,E,I1,I2,I3,R,D) u*I3;


% init
E   = 20;   
S   = N - E;
I1  = 0;
I2  = 0;
I3  = 0;
R   = 0;
D   = 0;

% capasity 
cap_ICU = inf;


% prpoergate
SL_positive = [1 2 3 6 11 19 29 42 53 66 72];

for itr=1:length(SL_positive)+7  % for 1000 days
    
    track(itr,:)  = [S E I1 I2 I3 R D];
    display(round([S E I1 I2 I3 R D sum([S E I1 I2 I3 R D])]))
        
    dS_now  = dS(S,E,I1,I2,I3,R,D); 
    dE_now  = dE(S,E,I1,I2,I3,R,D); 
    dI1_now = dI1(S,E,I1,I2,I3,R,D); 
    dI2_now = dI2(S,E,I1,I2,I3,R,D); 
    dI3_now = dI3(S,E,I1,I2,I3,R,D); 
    dR_now  = dR(S,E,I1,I2,I3,R,D); 
    dD_now  = dD(S,E,I1,I2,I3,R,D); 
    
    S       = S  + dS_now;
    E       = E  + dE_now;
    I1      = I1 + dI1_now;
    I2      = I2 + dI2_now;
    
    if I3<cap_ICU
        I3      = I3 + dI3_now;
        R       = R  + dR_now;
        D       = D  + dD_now;
    else        
        I3      = I3;
        R       = R  + dR_now;
        D       = D  + dD_now + dI3_now;
    end
    
end

%% log plot
figure
track_rnd = round(track);
semilogy(track_rnd(:,2),'+-b','LineWidth',1.5);hold on      % E
semilogy(track_rnd(:,3),'s-k','LineWidth',1.5);hold on      % mild
semilogy(track_rnd(:,4),'x-k','LineWidth',1.5);hold on      % severe
semilogy(track_rnd(:,5),'O-k','LineWidth',1.5);hold on      % critical
semilogy(track_rnd(:,6),'+-g','LineWidth',1.5);hold on      % R
semilogy(track_rnd(:,end),'O-r','LineWidth',1.5);hold on    % D
legend('Exposed','Mild','Severe','Critical','Recovered','Dead')
title('Model Predictions (with no Intervention)')
xlabel('Days');
ylabel('Number of Individuals');
set(gca,'fontsize',20);
saveas(gcf,'./Model_predictions_semiLog.tif');

figure;
semilogy(sum(track(:,3:end),2),'+-b','LineWidth',1.5);hold on   % mild + severe + critical
semilogy(SL_positive,'s-m','LineWidth',1.5);hold off            % SL+
legend('Model(All Infected)','Tested +ve in SL')
title('SL vs. Model (with no Intervention)')
xlabel('Days');
ylabel('Number of Individuals');
set(gca,'fontsize',20);
saveas(gcf,'./Model_vs_SL_semiLog.tif');



%% normla plot
figure
track_rnd = round(track);
plot(track_rnd(:,1),'--b','LineWidth',1.5);hold on      % S
plot(track_rnd(:,2),'+-b','LineWidth',1.5);hold on      % E
plot(track_rnd(:,3),'s-k','LineWidth',1.5);hold on      % mild
plot(track_rnd(:,4),'x-k','LineWidth',1.5);hold on      % severe
plot(track_rnd(:,5),'O-k','LineWidth',1.5);hold on      % critical
plot(track_rnd(:,6),'+-g','LineWidth',1.5);hold on      % R
plot(track_rnd(:,end),'O-r','LineWidth',1.5);hold on    % D
legend('Susceptible','Exposed','Mild','Severe','Critical','Recovered','Dead')
title('Model Predictions (with no Intervention)')
xlabel('Days');
ylabel('Number of Individuals');
set(gca,'fontsize',20);
saveas(gcf,'./Model_predictions.tif');

figure;
plot(sum(track(:,3:end),2),'+-b','LineWidth',1.5);hold on   % mild + severe + critical
plot(SL_positive,'s-m','LineWidth',1.5);hold off            % SL+
legend('Model(All Infected)','Tested +ve in SL')
title('SL vs. Model(with no Intervention)')
xlabel('Days');
ylabel('Number of Individuals');
set(gca,'fontsize',20);
saveas(gcf,'./Model_vs_SL.tif');




