

classdef START_NODE < handle
   properties 
        N               % [AU]      vector with an entry for each layer of the node, [unreported layer, reported layer, quarantined layer]
   end
         
   methods   
       function obj = START_NODE(N0)
           obj.N        = N0;       
       end
      
       function dN_out  = evolve(self,beta_effective,dt)                            % dt = updating rate [Days]    
           dN_out       = self.N * beta_effective * dt;                        % beta_effective [Day-1]  
           self.N       = self.N - dN_out;                                                  
       end       
   end
   
   
end