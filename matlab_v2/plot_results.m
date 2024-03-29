
function plot_results(evolutions, gt, t_evolve,dates)

FS = 14;

tvec = [1:t_evolve+1]';
pred_diagnosed = (evolutions.I0(2,:) + evolutions.I1(2,:) + evolutions.I2(2,:) + ...
                  evolutions.I3(2,:) + evolutions.R(2,:) + evolutions.D(2,:) )';

hold off
subplot(1,3,1);
plot(tvec, round(pred_diagnosed), 'k-x', gt.tvec, gt.confirmed, 'b-O')
set(gca,'fontsize',FS);
xlabel('Time [days]');
ylabel('No. of diagnosed cases');
xticks(tvec)
xticklabels(dates)
xtickangle(90)
legend('Predicted','Observed','location','NW');

subplot(1,3,2);
plot(tvec, round(sum(evolutions.D,1))', 'k-x', gt.tvec, gt.deaths, 'b-O')
set(gca,'fontsize',FS);
xlabel('Time [days]');
ylabel('No. of fatalities');
xticks(tvec)
xticklabels(dates)
xtickangle(90)
legend('Predicted','Observed','location','NW');


subplot(1,3,3);
plot(tvec, round(evolutions.R(2,:))', 'k-x', gt.tvec, gt.recovered, 'b-O')
set(gca,'fontsize',FS);
xlabel('Time [days]');
ylabel('No. of recoveries');
xticks(tvec)
xticklabels(dates)
xtickangle(90)
legend('Predicted','Observed','location','NW');


set(gcf, 'Units', 'Inches', 'Position', [1,1,14,6])

end