classdef (Abstract) Function1d

       methods
           function obj = Function1d(), end % constructor
       end % methdds
       

% Abstract methods used to define a common interface for all subclasses. 
% Abstract methods must be implemented by all subclasses.

   methods(Abstract = true, Static)
         v = F(x);     %  function definition
         d = x1(x);   %  1st derivative
         d = x2(x);   %  2nd derivative
         d = x3(x);  %  3rd derivative
         d = x4(x);  %  4th derivative
   end  % abstract methods
   
end % class