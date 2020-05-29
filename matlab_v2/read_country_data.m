

function [country_names, country_data] = read_country_data(fName)

    Data_string     = fileread(fName);
    ind_openSqBrace     = find(Data_string=='[');
    ind_closeSqBrace    = find(Data_string==']');

    % first country name
    next_country_name   = strtrim(Data_string(1:ind_openSqBrace(1)-1));
    ind_quotes          = find(next_country_name=='"');
    next_country_name   = next_country_name(ind_quotes(1)+1:ind_quotes(2)-1);

    for i=1:length(ind_openSqBrace)
        country_ind = i;
        country_names{i}    = next_country_name;
        country_dataStr     = strtrim(Data_string(ind_openSqBrace(country_ind):ind_closeSqBrace(country_ind)));
        country_data{i}     = jsondecode(country_dataStr);

        % next country name
        if ~(i==length(ind_openSqBrace))
            next_country_name   = strtrim(Data_string(ind_closeSqBrace(country_ind)+1:ind_openSqBrace(country_ind+1)-1));
            ind_quotes          = find(next_country_name=='"');
            next_country_name   = next_country_name(ind_quotes(1)+1:ind_quotes(2)-1);
        end
    end
    
end

