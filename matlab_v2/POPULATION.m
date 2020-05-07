
classdef POPULATION < handle
   properties 
       S        % 
       E0       %
       E1       %
       I0       % I0.N = [I0_unreported I0_reported]
       I1       %   "
       I2       %   "
       I3       %   "
       R        %   "
       D        %   "
       
       S_0      % ~ entire population
       dt       % evolution rate
       Beta_gov % governing transmission ratio vector
       Beta_eff % effective transmission ratio scaler  
       
       c_I0     % reporting rate c_I0(t) 
       c_I1     % reporting rate c_I1(t) 
       c_I2     % reporting rate c_I2(t) 
       c_I3     % reporting rate c_I3(t) 
   end
         
   methods   
       function obj = POPULATION(S_0,T,fr,Beta_gov,c,dt)                    % T.XX  = T_duration in the node XX
         obj.S          = S_NODE(S_0);                                  % fr.XX = fraction flowing in from the previous node(s)
         obj.E0         = INTERMEDIATE_NODE(0,      T.E0,    fr.E0);        % rr.XX = reporting_rate    
         obj.E1         = INTERMEDIATE_NODE(0,      T.E1,    fr.E1);
         obj.I0         = INTERMEDIATE_NODE([0 0],  T.I0,    fr.I0);
         obj.I1         = INTERMEDIATE_NODE([0 0],  T.I1,    fr.I1);
         obj.I2         = INTERMEDIATE_NODE([0 0],  T.I2,    fr.I2);
         obj.I3         = INTERMEDIATE_NODE([0 0],  T.I3,    fr.I3);
         obj.R          = R_NODE([0 0],[fr.R_I0 fr.R_I1 fr.R_I2 fr.R_I3]);
         obj.D          = D_NODE([0 0],fr.D);
         
         obj.S_0        = S_0;
         obj.dt         = dt;
         obj.Beta_gov   = Beta_gov;
         
         obj.c_I0       = c(:,1);
         obj.c_I1       = c(:,2);
         obj.c_I2       = c(:,3);
         obj.c_I3       = c(:,4);
       end
              
       function reset(self)
         self.S.N       = self.S_0;   
         self.E0.N      = 0;   
         self.E1.N      = 0;
         self.I0.N      = [0 0];
         self.I1.N      = [0 0];
         self.I2.N      = [0 0];
         self.I3.N      = [0 0];
         self.R.N       = [0 0];
         self.D.N       = [0 0];
       end
       
       function self = set_Beta_effective(self,day_idx)
         self.Beta_eff  = self.E1.N(1) * self.Beta_gov(day_idx) + ...
                          self.I0.N(1) * self.Beta_gov(day_idx) + ...
                          self.I1.N(1) * self.Beta_gov(day_idx) + ...                          
                          self.I2.N(1) * self.Beta_gov(day_idx) + ...
                          self.I3.N(1) * self.Beta_gov(day_idx);
                      
                      % self.I1.N(2) * self.Beta_gov(day_idx) + ...   % Edited!!
       end
       
       function pop_history = evolve(self,N_days)           
          for t = 1:N_days
              self.set_Beta_effective(t);
                                          
              for i = self.dt:self.dt:1
%                   self.I0.report(self.dt,self.c_I0(t));
%                   self.I1.report(self.dt,self.c_I1(t));
%                   self.I2.report(self.dt,self.c_I2(t));
%                   self.I3.report(self.dt,self.c_I3(t));

                  dN_out_S  = self.S.evolve(self.Beta_eff,self.dt);
                  dN_out_E0 = self.E0.evolve( dN_out_S      ,self.dt);
                  dN_out_E1 = self.E1.evolve( dN_out_E0     ,self.dt);
                  dN_out_I0 = self.I0.evolve([dN_out_E1 0]  ,self.dt);          % [in4I0_unreported in4I0_reported] = [dN_out_E1 0]
                  dN_out_I1 = self.I1.evolve([dN_out_E1 0]  ,self.dt);
                  dN_out_I2 = self.I2.evolve( dN_out_I1     ,self.dt);
                  dN_out_I3 = self.I3.evolve( dN_out_I2     ,self.dt);
                                    
                  self.R.evolve(dN_out_I0,...
                                dN_out_I1,...
                                dN_out_I2,...
                                dN_out_I3);
                  self.D.evolve(dN_out_I3);                                
              end
              self.I0.report(1,self.c_I0(t));
              self.I1.report(1,self.c_I1(t));
              self.I2.report(1,self.c_I2(t));
              self.I3.report(1,self.c_I3(t));

              
              pop_history.S(t)     = self.S.N;      
              pop_history.E0(t)    = self.E0.N;
              pop_history.E1(t)    = self.E1.N;
              pop_history.I0(:,t)  = self.I0.N;
              pop_history.I1(:,t)  = self.I1.N;
              pop_history.I2(:,t)  = self.I2.N;
              pop_history.I3(:,t)  = self.I3.N;
              pop_history.R(:,t)   = self.R.N;
              pop_history.D(:,t)   = self.D.N;
              
%               disp(round(cat(2,[pop_history.S(t);0],...
%                          [pop_history.E0(t);0],...
%                          [pop_history.E1(t);0],...
%                           pop_history.I0(:,t),...
%                           pop_history.I1(:,t),...
%                           pop_history.I2(:,t),...
%                           pop_history.I3(:,t),...
%                           pop_history.D(:,t),...
%                           pop_history.R(:,t))))
          end           
       end       
   end      
end











