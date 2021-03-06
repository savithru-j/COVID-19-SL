clear
clc
close all


num_threads = 8;
trueCan = false;
data = readOptData('srilanka',1:num_threads,trueCan);

% load tempcosts
[sortedCost inds_selected] = sort(data.cost);
figure(1);
plot(sortedCost)
xlabel('Realization');
ylabel('Cost');

max_sel = min(24, length(sortedCost));
inds_selected = inds_selected(1:max_sel);
sortedCost    = sortedCost(1:max_sel);

%inds_selected = (sum(data.pred_inf_unreported,1)<20e3 & sum(data.pred_inf_unreported,1)>3e3);
% inds_selected = find(sum(data.pred_inf_unreported,1)>6e3);
% inds_selected = find(sum(data.pred_inf_unreported,1)<inf);

data.pred_conf = data.pred_conf(:,inds_selected);
data.pred_recov = data.pred_recov(:,inds_selected);
data.pred_fatal = data.pred_fatal(:,inds_selected);
data.pred_inf_unreported = data.pred_inf_unreported(:,inds_selected);
data.pred_recov_unreported = data.pred_recov_unreported(:,inds_selected);
data.pred_fatal_unreported = data.pred_fatal_unreported(:,inds_selected);
data.beta = data.beta(:,inds_selected);
data.c0 = data.c0(:,inds_selected);
data.c1 = data.c1(:,inds_selected);
data.c2 = data.c2(:,inds_selected);
data.c3 = data.c3(:,inds_selected);
data.last = data.last(:,inds_selected);

plot_log = false;
plot_cumulative = true;

figure(2)
subplot(2,3,1)
hold off
if (plot_log)
    if (plot_cumulative)
        semilogy(data.obs_conf, 'k--'); hold on
        semilogy(data.pred_conf, '-');
    else
        semilogy(diff(data.obs_conf), 'k--'); hold on
        semilogy(diff(data.pred_conf), '-');
    end
else
    if (plot_cumulative)
        plot(data.obs_conf, 'k--'); hold on
        plot(data.pred_conf, '-');
    else
        plot(diff(data.obs_conf), 'k--'); hold on
        plot(diff(data.pred_conf), '-');
    end
end
xlabel('Day')
ylabel('No. of confirmed cases')

subplot(2,3,2)
hold off
if (plot_log)
    semilogy(data.obs_recov, 'k--'); hold on
    semilogy(data.pred_recov, '-')
else
    plot(data.obs_recov, 'k--'); hold on
    plot(data.pred_recov, '-')
end
xlabel('Day')
ylabel('No. of recovered cases')
% xlim([0,65])

subplot(2,3,3)
hold off
if (plot_log)
    semilogy(data.obs_fatal, 'k--'); hold on
    semilogy(data.pred_fatal, '-')
else
    plot(data.obs_fatal, 'k--'); hold on
    plot(data.pred_fatal, '-')
end
xlabel('Day')
ylabel('No. of fatalities')

subplot(2,3,4)
hold off
if (plot_log)
    semilogy(data.pred_inf_unreported, '-')
else
    plot(data.pred_inf_unreported, '-')
end
xlabel('Day')
ylabel('No. of infected-unreported cases')

subplot(2,3,5)
hold off
if (plot_log)
    semilogy(data.pred_recov_unreported, '-')
else
    plot(data.pred_recov_unreported, '-')
end
xlabel('Day')
ylabel('No. of recovered-unreported cases')
% xlim([0,65])

subplot(2,3,6)
hold off
if (plot_log)
    semilogy(data.pred_fatal_unreported, '-')
else
    plot(data.pred_fatal_unreported, '-')
end

xlabel('Day')
ylabel('No. of unreported-fatalities')

set(gcf, 'Units', 'Inches', 'Position', [1,3,16,5])

%%
figure(3)
subplot(2,2,1)
denom = max(sortedCost) - min(sortedCost);
denom(denom == 0) = 1;
cost_color = (sortedCost - min(sortedCost)) ./ denom;
plotSorted(data.beta,cost_color);
if trueCan
  hold on; plot(data.beta_true,'r--','linewidth',2);
end
% plot(mean(data.beta,2),'b--','linewidth',2);hold off
xlabel('Day')
ylabel('Beta')
% legend('Optimal1','Optimal2','Optimal3','Optimal4','Optimal5')

subplot(2,2,2)
plotSorted(data.c0,cost_color);
if trueCan
  hold on; plot(data.c0_true,'r--','linewidth',2); hold off
end
xlabel('Day')
ylabel('c_0 = c_e')

subplot(2,2,3)
plotSorted(data.c1,cost_color);
if trueCan
   hold on; plot(data.c1_true,'r--','linewidth',2); hold off
end
xlabel('Day')
ylabel('c_1')

subplot(2,2,4)
plotSorted(data.c2,cost_color);
if trueCan
  hold on; plot(data.c2_true,'r--','linewidth',2); hold off
end
xlabel('Day')
ylabel('c_2')

% subplot(3,2,5)
% plotSorted(data.c3,cost_color);
% xlabel('Day')
% ylabel('c_3')

set(gcf, 'Units', 'Inches', 'Position', [1,3,14,8])


function data = readOptData(country, seed_range, trueCan)

if trueCan
  true_params_file = sprintf('csv_data/%s_params.txt', country);
  data_true_params = importdata(true_params_file);
  nt = (size(data_true_params,1) - 10) / 5;
  data.beta_true = data_true_params(1:nt, :);
  data.c0_true = data_true_params((  nt+1):2*nt, :);
  data.c1_true = data_true_params((2*nt+1):3*nt, :);
  data.c2_true = data_true_params((3*nt+1):4*nt, :);
  data.c3_true = data_true_params((4*nt+1):5*nt, :);
  data.last_true = data_true_params(end-10:end,:);
end

observed_file = sprintf('csv_data/%s.txt', country);
data_obs = importdata(observed_file);
data_obs = reshape(data_obs(2:end),3,[])';
data.obs_conf = data_obs(:,1);
data.obs_recov = data_obs(:,2);
data.obs_fatal = data_obs(:,3);


data_params = [];
data_pred = [];

for seed = seed_range

    params_file = sprintf('../C++/build/release/results/%s_params_seed%d.txt', country, seed);
    pred_file = sprintf('../C++/build/release/results/%s_prediction_seed%d.txt', country, seed);

    data_params = [data_params, importdata(params_file)];
    data_pred = [data_pred, importdata(pred_file)];
    
end

data.pred_conf = data_pred(:,1:6:end);
data.pred_recov = data_pred(:,2:6:end);
data.pred_fatal = data_pred(:,3:6:end);
data.pred_inf_unreported = data_pred(:,4:6:end);
data.pred_recov_unreported = data_pred(:,5:6:end);
data.pred_fatal_unreported = data_pred(:,6:6:end);

nt = (size(data_params,1) - 10 - 1) / 5
data.beta = data_params(1:nt, :);
data.c0 = data_params((  nt+1):2*nt, :);
data.c1 = data_params((2*nt+1):3*nt, :);
data.c2 = data_params((3*nt+1):4*nt, :);
data.c3 = data_params((4*nt+1):5*nt, :);
data.last = data_params(end-10-1:end-1,:);
data.cost = data_params(end,:);

end


function plotSorted(pltData,cost_color)
    hold on
    for i=length(cost_color):-1:1
        plot(pltData(:,i), '-', 'Color', [cost_color(i) cost_color(i) cost_color(i)])
    end
    hold off
end


