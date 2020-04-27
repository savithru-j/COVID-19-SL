
classdef D_NODE < handle
   properties 
        N               % [AU]      vector with an entry for each layer of the node, [unreported layer, reported layer, quarantined layer]
        frac_in         % [AU]      fractions flowing in from the previous node        
   end
         
   methods   
       function obj = D_NODE(N0,fraction_from_previous_node)
           obj.N        = N0;       
           obj.frac_in  = fraction_from_previous_node(1);                     % from I3 
       end
      
       function evolve(self,dN_out_from_previous_node)        % dt = updating rate [Days]               
           self.N       =   self.N ...
                          + dN_out_from_previous_node * self.frac_in; 
       end       
   end
   
   
end