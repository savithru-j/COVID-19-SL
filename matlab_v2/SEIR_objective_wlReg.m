
function [f evolutions]  = SEIR_objective_wlReg(x,model_population,gt,err_type)
    
    %gt.confirmed,deaths,recovered
    
    t_evolve                    = length(gt.confirmed);
    x                           = reshape(x,t_evolve,5);
    model_population.Beta_gov   = x(:,1)/model_population.S.N;
    model_population.c_I0       = x(:,2);
    model_population.c_I1       = x(:,3);
    model_population.c_I2       = x(:,4);
    model_population.c_I3       = x(:,5);
        
    model_population.reset;
    model_population.E0.N = 5;
    evolutions = model_population.evolve(t_evolve);
    
%     plot_results(evolutions, gt, t_evolve,[1:t_evolve]);
%     drawnow
    
    model.confirmed = evolutions.I0(2,:)'+evolutions.I1(2,:)'+evolutions.I2(2,:)'+evolutions.I3(2,:)'+evolutions.R(2,:)'+evolutions.D(2,:)';
    model.deaths    = evolutions.D(2,:)';
    model.recovered = evolutions.R(2,:)';

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
    
     %% 


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
    
    fprintf('conf = %1.4e | dead = %1.4e | rec = %1.4e | regs = [%1.4e %1.4e %1.4e %1.4e %1.4e] | cost = %1.4e\n', ...
            w(1)*f1, w(2)*f2, w(3)*f3,  w(4)*reg_beta, w(5)*reg_c_I0, w(6)*reg_c_I1, w(7)*reg_c_I2, w(8)*reg_c_I3, f)        
end



