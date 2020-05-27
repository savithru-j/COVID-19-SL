clear
clc
close all

country = 'srilanka64';
suffix = '';

observed_file = sprintf('csv_data/%s.txt', country);
params_file = sprintf('../C++/build/release/results/%s_params%s.txt', country, suffix);
pred_file = sprintf('../C++/build/release/results/%s_prediction%s.txt', country, suffix);

data_obs = importdata(observed_file);
data_obs = reshape(data_obs(2:end),3,[])';

data_params = importdata(params_file);
data_pred = importdata(pred_file);

pred_conf = data_pred(:,1:3:end);
pred_recov = data_pred(:,2:3:end);
pred_fatal = data_pred(:,3:3:end);

nt = (size(data_params,1) - 10) / 3
beta = data_params(1:nt, :);
ce = data_params((  nt+1):2*nt, :);
c1 = data_params((2*nt+1):3*nt, :);

data_params(end-9:end,:)

figure(1)
subplot(1,3,1)
hold on
plot(data_obs(:,1), 'k--')
plot(pred_conf, '-')
hold off
xlabel('Day')
ylabel('No. of confirmed cases')

subplot(1,3,2)
hold on
plot(data_obs(:,2), 'k--')
plot(pred_recov, '-')
hold off
xlabel('Day')
ylabel('No. of recovered cases')

subplot(1,3,3)
hold on
plot(data_obs(:,3), 'k--')
plot(pred_fatal, '-')
hold off
xlabel('Day')
ylabel('No. of fatalities')

set(gcf, 'Units', 'Inches', 'Position', [1,3,16,5])

%%
figure(2)
subplot(3,1,1)
plot(beta, '-')
xlabel('Day')
ylabel('Beta')
legend('Optimal1','Optimal2','Optimal3','Optimal4','Optimal5')

subplot(3,1,2)
plot(ce, '-')
xlabel('Day')
ylabel('c_e')

subplot(3,1,3)
plot(c1, '-')
xlabel('Day')
ylabel('c_1')

set(gcf, 'Units', 'Inches', 'Position', [1,3,10,8])

