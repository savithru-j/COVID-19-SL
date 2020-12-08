clear
clc

A = importdata('SL_district_data.csv');

f = fopen('SL_district_data.js','w');

fprintf(f, 'let SL_district_data = \n{\n');

for i=1:26
   fprintf(f, '\t"%s": [\n', A.textdata{1,i+1});
   
   for j=1:(size(A.textdata,1)-1)
        fprintf(f, '\t\t{x: "%s", y: %d},\n', A.textdata{j+1,1}, A.data(j,i));
   end
   fprintf(f, '\t],\n');
end

fprintf(f,'}');
fclose(f);