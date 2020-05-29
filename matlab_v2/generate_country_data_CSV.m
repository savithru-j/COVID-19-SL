clear
clc
close all

[country_names, country_data] = read_country_data('../world_data.js');

folder_path = "csv_data/";

countries_to_export = ["Sri Lanka", "US"];
populations = [21323733, 329064917];

for i = 1:length(countries_to_export)
    
    data = country_data{country_names == countries_to_export(i)};
    
    filename = folder_path + lower(erase(countries_to_export(i)," ")) + ".txt";
    f = fopen(filename, 'w');
    fprintf(f, '%d\n', populations(i));
    
    for j = 1:length(data)
       if (data(j).confirmed > 0)
           fprintf(f, '%d %d %d\n', data(j).confirmed, data(j).recovered, data(j).deaths);
       end
    end
    
    fclose(f);
end