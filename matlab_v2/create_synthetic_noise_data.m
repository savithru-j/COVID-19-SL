clear
clc

rng(1);

for i = 1:1
   data = importdata(sprintf('csv_data/synthetic/synthetic%d.txt',i));
   N = data(1);
   data = reshape(data(2:end),3,[])';
   
%    plot(data(:,1));
%    hold on
   
   for j = 1:50
        obs_data = poissrnd(data);
%         plot(obs_data(:,1));

        f = fopen(sprintf('csv_data/synthetic%d_noise%d.txt',i,j), 'w');
        fprintf(f,'%d\n', N);
        fprintf(f,'%d %d %d\n',obs_data');
        fclose(f);
   end
end

