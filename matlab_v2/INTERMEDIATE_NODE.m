
classdef INTERMEDIATE_NODE < handle
   properties 
        N               % [AU]      vector with an entry for each layer of the node, [unreported layer, reported layer, quarantined layer]
        T               % [Days]    time duration spent in the node
        frac_in         % [AU]      fraction flowing in from the previous node
        reporting_rate  % [Day-1]   rate at which the unreported cases gets reorted for this node
   end
         
   methods   
       function obj = INTERMEDIATE_NODE(N0,T_duration,fraction_from_previous_node,reporting_rate)
           obj.N                = N0;
           obj.T                = T_duration;
           obj.frac_in          = fraction_from_previous_node;
           obj.reporting_rate   = reporting_rate;           
       end
      
       function dN_out = evolve(self,dN_out_from_previous_node,dt)                            % dt = updating rate [Days]                          
%            if length(self.N)>1
%                self.report(dt);
%            end           

           dN_out       =   self.N * (1/self.T) * dt;
           self.N       =   self.N ...
                          + dN_out_from_previous_node * self.frac_in ...
                          - dN_out;                                             
       end
       
       function report(self,dt)
           dN_reported  = self.N(1) * self.reporting_rate * dt;
           self.N(1)    = self.N(1) - dN_reported;
           self.N(2)    = self.N(2) + dN_reported;
       end
   end
   
   
end