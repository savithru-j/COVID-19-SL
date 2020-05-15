
function [f evolutions]  = SEIR_objective_wlReg(x,model_population,gt,err_type)
        
    t_evolve                    = length(gt.confirmed);                     % gt struture" gt.confirmed, gt.deaths, gt.recovered
    
    %% set trainables from x
    x                               = reshape(x,t_evolve,6);
    model_population.Beta_gov       = x(:,1)/model_population.S.N;
    model_population.c_I0           = x(:,2);
    model_population.c_I1           = x(:,3);
    model_population.c_I2           = x(:,4);
    model_population.c_I3           = x(:,5);
    
    model_population.E0.T           = x(1,6);
    model_population.E1.T           = x(2,6);
    model_population.I0.T           = x(3,6);
    model_population.I1.T           = x(4,6);
    model_population.I2.T           = x(5,6);
    model_population.I3.T           = x(6,6);

    model_population.I0.frac_in     = x(7,6);                                                               % fr.I0   = 0.3;
    model_population.I1.frac_in     = 1 - model_population.I0.frac_in;                                      % fr.I1   = 1 - fr.I0;
    model_population.R.frac_in_I1   = x(8,6);                                                               % fr.R_I1 = 0.8;                 % 80% of I1 cases recover
    model_population.I2.frac_in     = 1 - model_population.R.frac_in_I1;                                    % fr.I2   = 1 - fr.R_I1;         % rest goes to I2    
    model_population.R.frac_in_I2   = x(9,6); %/model_population.I2.frac_in;                                   % fr.R_I2 = 0.15/fr.I2;          % 15% of sever cases recover
    model_population.I3.frac_in     = 1 - model_population.R.frac_in_I2;                                    % fr.I3   = 1 - fr.R_I2;
    model_population.D.frac_in      = min(x(10,6)/(model_population.I2.frac_in*model_population.I3.frac_in), 1);    % fr.D    = 0.02/(fr.I2*fr.I3);  % 2% fatalify rate    
    model_population.R.frac_in_I3   = 1 - model_population.D.frac_in;                                       % fr.R_I3 = 1 - fr.D;   
    
    %Alternate way:
    %model_population.R.frac_in_I3   = x(10,6);
    %model_population.D.frac_in      = 1 - model_population.R.frac_in_I3;    
    
%     if (model_population.R.frac_in_I3 < 0)
%        disp('negative frac');     
%     end
    
    model_population.reset;
    
    model_population.E0.N = 5;
    model_population.R.N(2) = 1; %Special case for Sri Lanka
    model_population.S.N = model_population.S_0 - model_population.E0.N - model_population.R.N(2);
    
    evolutions = model_population.evolve(t_evolve);    
%     plot_results(evolutions, gt, t_evolve,[1:t_evolve]);                  % plot evaluation  
%     drawnow
    
    ind = 2:size(evolutions.I0, 2); %Skip initial condition
    model.confirmed = (evolutions.I0(2,ind) + evolutions.I1(2,ind) + evolutions.I2(2,ind) + ...
                       evolutions.I3(2,ind) + evolutions.R(2,ind) + evolutions.D(2,ind))';
    model.deaths    = evolutions.D(2,ind)';
    model.recovered = evolutions.R(2,ind)';

    switch err_type
        case 'log_type1'                
            delta = 1e-5;
            err_c = log(model.confirmed+delta) - log(gt.confirmed'+delta);               % log error type1 (doesn't work great)
            err_d = log(model.deaths+delta)    - log(gt.deaths'+delta);
            err_r = log(model.recovered+delta) - log(gt.recovered'+delta);
            f1  = (err_c'*err_c) / (log(gt.confirmed'+delta)'*log(gt.confirmed'+delta));  % L2 sq. confirmed
            f2  = (err_d'*err_d) / (log(gt.deaths'+delta)'*log(gt.deaths'+delta));        % L2 sq. deaths
            f3  = (err_r'*err_r) / (log(gt.recovered'+delta)'*log(gt.recovered'+delta));  % L2 sq. recovered   
        case 'log_type2'
            delta = 1e-5;
            err_c = log(abs(model.confirmed - gt.confirmed'+delta));                      % log error type2 (doesn't work great)
            err_d = log(abs(model.deaths    - gt.deaths'+delta));
            err_r = log(abs(model.recovered - gt.recovered'+delta));            
            f1  = (err_c'*err_c) / (log(gt.confirmed'+delta)'*log(gt.confirmed'+delta));  % L2 sq. confirmed
            f2  = (err_d'*err_d) / (log(gt.deaths'+delta)'*log(gt.deaths'+delta));        % L2 sq. deaths
            f3  = (err_r'*err_r) / (log(gt.recovered'+delta)'*log(gt.recovered'+delta));  % L2 sq. recovered               
        case 'L2_type1'
            err_c = model.confirmed - gt.confirmed';              % normal L2 error (works ok)               
            err_d = model.deaths    - gt.deaths';
            err_r = model.recovered - gt.recovered';
            f1  = (err_c'*err_c) / (gt.confirmed*gt.confirmed');  % L2 sq. confirmed
            f2  = (err_d'*err_d) / (gt.deaths*gt.deaths');        % L2 sq. deaths
            f3  = (err_r'*err_r) / (gt.recovered*gt.recovered');  % L2 sq. recovered   
        case 'L2_type2'                  
            err_c = diff(model.confirmed) - diff(gt.confirmed)';  % MSE like L2 error (generates NaNs!)  
            err_d = diff(model.deaths)    - diff(gt.deaths)';
            err_r = diff(model.recovered) - diff(gt.recovered)';
            f1  = (err_c'*err_c) / t_evolve;  % L2 sq. confirmed
            f2  = (err_d'*err_d) / t_evolve;  % L2 sq. deaths
            f3  = (err_r'*err_r) / t_evolve;  % L2 sq. recovered   
    end

    %% wavelet regularizer
    n_wl    = fix(log2(t_evolve));    
    [c l]   = wavedec(x(:,1),n_wl,'db1');% 'haar','db1','sym2'
    reg_beta= norm(c,1);
    [c l]   = wavedec(x(:,2),n_wl,'db1');% 'haar','db1','sym2'
    reg_c_I0= norm(c,1);
    [c l]   = wavedec(x(:,3),n_wl,'db1');% 'haar','db1','sym2'
    reg_c_I1= norm(c,1);
    [c l]   = wavedec(x(:,4),n_wl,'db1');% 'haar','db1','sym2'
    reg_c_I2= norm(c,1);
    [c l]   = wavedec(x(:,5),n_wl,'db1');% 'haar','db1','sym2'
    reg_c_I3= norm(c,1);
    
    w_reg = 0.001;
%              L2-conf    L2-dead     L2-rec      Reg_beta        Reg_c_I0 ...    
    w       = [1          2           0           w_reg           w_reg           w_reg           w_reg           w_reg];
    f       =  w(1)*f1  + w(2)*f2   + w(3)*f3   + w(4)*reg_beta + w(5)*reg_c_I0 + w(6)*reg_c_I1 + w(7)*reg_c_I2 + w(8)*reg_c_I3;
    
    if (isnan(f))
       disp('found NaN'); 
    end
    
    fprintf('conf = %1.4e | dead = %1.4e | rec = %1.4e | regs = [%1.4e %1.4e %1.4e %1.4e %1.4e] | cost = %1.4e\n', ...
            w(1)*f1, w(2)*f2, w(3)*f3,  w(4)*reg_beta, w(5)*reg_c_I0, w(6)*reg_c_I1, w(7)*reg_c_I2, w(8)*reg_c_I3, f)        
end



