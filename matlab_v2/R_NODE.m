
classdef R_NODE < handle
   properties 
        N               % [AU]      vector with an entry for each layer of the node, [unreported layer, reported layer, quarantined layer]
        frac_in_I0      % [AU]      fractions flowing in from the previous nodes (R has multiple input nodes)
        frac_in_I1
        frac_in_I2
        frac_in_I3
   end
         
   methods   
       function obj = R_NODE(N0,fraction_from_previous_nodes)
           obj.N            = N0;       
           obj.frac_in_I0   = fraction_from_previous_nodes(1);                     % order [I0 I1 I2 I3]
           obj.frac_in_I1   = fraction_from_previous_nodes(2);                     % order [I0 I1 I2 I3]
           obj.frac_in_I2   = fraction_from_previous_nodes(3);                     % order [I0 I1 I2 I3]
           obj.frac_in_I3   = fraction_from_previous_nodes(4);                     % order [I0 I1 I2 I3]
       end
      
       function evolve(self,dN_out_I0,dN_out_I1,dN_out_I2,dN_out_I3)                                       % dt = updating rate [Days]               
           self.N       =   self.N ...
                          + dN_out_I0 * self.frac_in_I0 ...
                          + dN_out_I1 * self.frac_in_I1 ...
                          + dN_out_I2 * self.frac_in_I2 ...
                          + dN_out_I3 * self.frac_in_I3 ;
  
       end       
   end
   
   
end