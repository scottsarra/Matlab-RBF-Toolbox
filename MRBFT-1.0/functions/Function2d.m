classdef (Abstract) Function2d

       methods
           function obj = Function2d(), end % constructor
       end % methdds
       

% Abstract methods used to define a common interface for all subclasses. 
% Abstract methods must be implemented by all subclasses.

   methods(Abstract = true, Static)
         v = F(x,y);     %  function definition
         d = x1(x,y);    %  1st derivative wrt x
         d = x2(x,y);    %  2nd derivative
         d = x3(x,y);    %  3rd derivative
         d = x4(x,y);    %  4th derivative
         d = y1(x,y);    %  1st derivative wrt y
         d = y2(x,y);    %  2nd derivative
         d = y3(x,y);    %  3rd derivative
         d = y4(x,y);    %  4th derivative
         d = G(x,y);     %  Gradient
         d = L(x,y);     %  Laplacian
         d = B(x,y);     %  Biharmonic
         d = p12(x,y);   %  mixed partials
         d = p21(x,y);
         d = p22(x,y)
   end  % abstract methods
   
end % class