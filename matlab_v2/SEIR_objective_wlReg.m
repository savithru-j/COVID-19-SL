
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

    f1  = norm(model.confirmed - gt.confirmed,2) / norm(gt.confirmed,2);        % L2 confirmed
    f2  = norm(model.deaths    - gt.deaths,2)    / norm(gt.deaths,   2);        % L2 deaths
    f3  = norm(model.recovered - gt.recovered,2) / norm(gt.recovered,2);        % L2 recovered   
    
    % wavelet regularizer
    n_wl    = fix(log2(t_evolve));    
    [c l]   = wavedec(x(1:t_evolve),n_wl,'haar');
    reg     = norm(c,1);

%              L2-conf    L2-dead     L2-rec      Reg  
    w       = [1          1           0           0.01];
    f       =  w(1)*f1  + w(2)*f2   + w(3)*f3   + w(4)*reg;
    
    disp(sprintf('conf = %d | dead = %d | rec = %d | reg = %d', w(1)*f1, w(2)*f2, w(3)*f3, w(4)*reg))        
end



