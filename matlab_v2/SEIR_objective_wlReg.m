
function [f evolutions]  = SEIR_objective_wlReg(x,model_population,gt)
    
    %gt.confirmed,deaths,recovered
    
    t_evolve                            = length(gt.confirmed);
    model_population.Beta_gov           = x(1:t_evolve)/model_population.S.N;
    model_population.I0.reporting_rate  = x(t_evolve+1);
    model_population.I1.reporting_rate  = x(t_evolve+2);
    model_population.I2.reporting_rate  = x(t_evolve+3);
    model_population.I3.reporting_rate  = x(t_evolve+4);
        
    model_population.reset;
    model_population.E0.N = 5;
    evolutions = model_population.evolve(t_evolve);
    
    model.confirmed = evolutions.I0(2,:)'+evolutions.I1(2,:)'+evolutions.I2(2,:)'+evolutions.I3(2,:)'+evolutions.R(2,:)'+evolutions.D(2,:)';
    model.deaths    = evolutions.D(2,:)';
    model.recovered = evolutions.R(2,:)';
    
    err_c = model.confirmed - gt.confirmed';
    err_d = model.deaths    - gt.deaths';
    err_r = model.recovered - gt.recovered';

    f1  = (err_c'*err_c) / (gt.confirmed*gt.confirmed');  % L2 sq. confirmed
    f2  = (err_d'*err_d) / (gt.deaths*gt.deaths');        % L2 sq. deaths
    f3  = (err_r'*err_r) / (gt.recovered*gt.recovered');  % L2 sq. recovered   
    
    % wavelet regularizer
    n_wl    = fix(log2(t_evolve));    
    [c l]   = wavedec(x(1:t_evolve),n_wl,'haar');
    reg     = norm(c,1);

%              L2-conf    L2-dead     L2-rec      Reg  
    w       = [1          1           0           0.01];
    f       =  w(1)*f1  + w(2)*f2   + w(3)*f3   + w(4)*reg;
    
    fprintf('conf = %1.4e | dead = %1.4e | rec = %1.4e | reg = %1.4e | cost = %1.4e\n', ...
            w(1)*f1, w(2)*f2, w(3)*f3, w(4)*reg, f)        
end



