% clear
% clc
% close all

num_threads = 40;
data = readOptData('synthetic_new1',1:num_threads);

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

col_ind = [1:5];

figure(11)
subplot(1,3,1)
hold off
semilogy(data.pred_conf(:,col_ind), '-')
hold on
semilogy(data.obs_conf, 'k--')
xlabel('Day')
ylabel('No. of confirmed cases')

subplot(1,3,2)
hold off
semilogy(data.pred_recov(:,col_ind), '-')
hold on
semilogy(data.obs_recov, 'k--')
xlabel('Day')
ylabel('No. of recovered cases')
% xlim([0,65])

subplot(1,3,3)
hold off
semilogy(data.pred_fatal(:,col_ind), '-')
hold on
semilogy(data.obs_fatal, 'k--')
xlabel('Day')
ylabel('No. of fatalities')

set(gcf, 'Units', 'Inches', 'Position', [1,3,16,5])

%%

figure(12)
subplot(3,2,1)
hold off
plot(data.beta(:,col_ind), '-')
hold on
plot(data.beta_true,'k--')
xlabel('Day')
ylabel('Beta')
% legend('Optimal1','Optimal2','Optimal3','Optimal4','Optimal5','True')
legend('Optimal1','True')

subplot(3,2,2)
hold off
plot(data.c0(:,col_ind), '-')
hold on
plot(data.c0_true,'k--')
xlabel('Day')
ylabel('c_0 = c_e')

subplot(3,2,3)
plot(data.c1(:,col_ind), '-')
hold on
plot(data.c1_true,'k--')
xlabel('Day')
ylabel('c_1')

subplot(3,2,4)
plot(data.c2(:,col_ind), '-')
hold on
plot(data.c2_true,'k--')
xlabel('Day')
ylabel('c_2')

subplot(3,2,5)
plot(data.c3(:,col_ind), '-')
hold on
plot(data.c3_true,'k--')
xlabel('Day')
ylabel('c_3')

set(gcf, 'Units', 'Inches', 'Position', [1,3,14,8])




function data = readOptData(country, seed_range)

observed_file = sprintf('csv_data/%s.txt', country);
true_params_file = sprintf('csv_data/%s_params.txt', country);
params_file = sprintf('../C++/build/release/results/%s_params_seed%s.txt', country, suffix);
pred_file = sprintf('../C++/build/release/results/%s_prediction_seed%s.txt', country, suffix);

data_obs = importdata(observed_file);
data_obs = reshape(data_obs(2:end),3,[])';
data.obs_conf = data_obs(:,1);
data.obs_recov = data_obs(:,2);
data.obs_fatal = data_obs(:,3);

data_true_params = importdata(true_params_file);
data_params = importdata(params_file);
data_pred = importdata(pred_file);

data.pred_conf = data_pred(:,1:6:end);
data.pred_recov = data_pred(:,2:6:end);
data.pred_fatal = data_pred(:,3:6:end);

nt = (size(data_params,1) - 11 -1) / 5;
data.beta = data_params(1:nt, :);
data.c0 = data_params((  nt+1):2*nt, :);
data.c1 = data_params((2*nt+1):3*nt, :);
data.c2 = data_params((3*nt+1):4*nt, :);
data.c3 = data_params((4*nt+1):5*nt, :);
data.last = data_params(end-10:end,:);

nt = (size(data_true_params,1) - 11) / 5;
data.beta_true = data_true_params(1:nt, :);
data.c0_true = data_true_params((  nt+1):2*nt, :);
data.c1_true = data_true_params((2*nt+1):3*nt, :);
data.c2_true = data_true_params((3*nt+1):4*nt, :);
data.c3_true = data_true_params((4*nt+1):5*nt, :);
data.last_true = data_true_params(end-10:end,:);

end