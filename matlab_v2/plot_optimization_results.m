clear
clc
close all

data = readOptData('srilanka64','');

% data0 = readOptData('srilanka64_T1','');
% data = readOptData('srilanka64_T14','');
% 
% data0.last
% data.last
% 
% data.pred_conf(:,1) = data0.pred_conf(:,1);
% data.pred_recov(:,1) = data0.pred_recov(:,1);
% data.pred_fatal(:,1) = data0.pred_fatal(:,1);
% data.beta(:,1) = data0.beta(:,1);
% data.c0(:,1) = data0.c0(:,1);
% data.c1(:,1) = data0.c1(:,1);
% data.c2(:,1) = data0.c2(:,1);
% data.c3(:,1) = data0.c3(:,1);

figure(1)
subplot(2,3,1)
hold on
plot(data.obs_conf, 'k--')
plot(data.pred_conf, '-')
plot(data.pred_inf_unreported + data.pred_recov_unreported + data.pred_fatal_unreported, '-.')
hold off
xlabel('Day')
ylabel('No. of confirmed cases')

subplot(2,3,2)
hold on
plot(data.obs_recov, 'k--')
plot([zeros(0,5); data.pred_recov], '-')
hold off
xlabel('Day')
ylabel('No. of recovered cases')
% xlim([0,65])

subplot(2,3,3)
hold on
plot(data.obs_fatal, 'k--')
plot(data.pred_fatal, '-')
hold off
xlabel('Day')
ylabel('No. of fatalities')

subplot(2,3,4)
plot(data.pred_inf_unreported, '-')
hold off
xlabel('Day')
ylabel('No. of confirmed cases')

subplot(2,3,5)
plot(data.pred_recov_unreported, '-')
xlabel('Day')
ylabel('No. of recovered cases')
% xlim([0,65])

subplot(2,3,6)
plot(data.pred_fatal_unreported, '-')

xlabel('Day')
ylabel('No. of fatalities')

set(gcf, 'Units', 'Inches', 'Position', [1,3,16,5])

%%
figure(2)
subplot(3,2,1)
plot(data.beta, '-')
xlabel('Day')
ylabel('Beta')
legend('Optimal1','Optimal2','Optimal3','Optimal4','Optimal5')

subplot(3,2,2)
plot(data.c0, '-')
xlabel('Day')
ylabel('c_0 = c_e')

subplot(3,2,3)
plot(data.c1, '-')
xlabel('Day')
ylabel('c_1')

subplot(3,2,4)
plot(data.c2, '-')
xlabel('Day')
ylabel('c_2')

subplot(3,2,5)
plot(data.c3, '-')
xlabel('Day')
ylabel('c_3')

set(gcf, 'Units', 'Inches', 'Position', [1,3,14,8])


function data = readOptData(country, suffix)

observed_file = sprintf('csv_data/%s.txt', country);
params_file = sprintf('../C++/build/release/results/%s_params%s.txt', country, suffix);
pred_file = sprintf('../C++/build/release/results/%s_prediction%s.txt', country, suffix);

data_obs = importdata(observed_file);
data_obs = reshape(data_obs(2:end),3,[])';
data.obs_conf = data_obs(:,1);
data.obs_recov = data_obs(:,2);
data.obs_fatal = data_obs(:,3);

data_params = importdata(params_file);
data_pred = importdata(pred_file);

data.pred_conf = data_pred(:,1:6:end);
data.pred_recov = data_pred(:,2:6:end);
data.pred_fatal = data_pred(:,3:6:end);
data.pred_inf_unreported = data_pred(:,4:6:end);
data.pred_recov_unreported = data_pred(:,5:6:end);
data.pred_fatal_unreported = data_pred(:,6:6:end);

nt = (size(data_params,1) - 11) / 5
data.beta = data_params(1:nt, :);
data.c0 = data_params((  nt+1):2*nt, :);
data.c1 = data_params((2*nt+1):3*nt, :);
data.c2 = data_params((3*nt+1):4*nt, :);
data.c3 = data_params((4*nt+1):5*nt, :);
data.last = data_params(end-10:end,:);

end