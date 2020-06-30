clear
clc
close all

rng(1)

nt = 100;

last_params = [3; 2; 6; 6; 4; 10; 0.3; 0.8; 0.75; 0.02; 14];

param_vecs = [];

Npop = 1e7; %population

for i = 1:10
   betaN = getrandomcurve(nt, 0.05, 0.95);
   c0 = zeros(size(betaN));
   c1 = getrandomcurve(nt, 0.05, 0.95);
   c2 = ones(size(betaN));
   c3 = c2;
   vec = [betaN; c0; c1; c2; c3; last_params];
   param_vecs = [param_vecs, vec];
   
   writematrix(vec, sprintf('csv_data/synthetic%d_params.txt',i),'Delimiter',' ')
end

plot(param_vecs)
% plot(param_vecs(nt+1:2*nt,:))
% ylim([0,2])

function y = getrandomcurve(nt, ymin, ymax)

nc = round(0.05*nt);

while (true)
    mean = ymin + (ymax-ymin)*rand();
    c = zeros(nt,1);
    ind = randi(round(0.10*nt), nc, 1);
    c(ind,:) = rand(nc,1);
    c(1,:) = mean*sqrt(nt);
    y = idct(c);
    if (ymin <= min(y) && max(y) <= ymax)
        break;
    end
end

end