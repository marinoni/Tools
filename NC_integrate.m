function [z]=NC_integrate(varargin)

%This function performs an integral on the input vector
%with the Newton-Cotes closed formulas of order 1,2,3 or 4
%
%By using higher orders, the grid step is accordingly increased. 
%E.g., 1st order takes as grid step the actual step.
%      2nd order takes as step the sum of two adiacent grid steps.
%Therefore, depending on the original grid step, higher order methods
%may not always be more accurate lower order methods.
%In particular, this is usually the case for order three.
%
%Formulas for higher order terms are derived only in case of equally
%spaced grids in each column. If the grid(s) are not equally spaced, the
%first order method should be used, as its formulation does not depend
%on the grid spacing. This will be upgraded in future versions
%
%INPUT:
%
%-y: vector to be integrated. For matrices, the integral is perfomed
%    column-wise
%-x: radial vector where y is defined. Default: equispaced vector
%-n: order of the algorithm. Default: 1 
%    1: Trapezoid rule
%    2: Cavalieri-Simpson rule
%    3: Simpson's 3/8 rule
%    4: Boole's rule
%-s: if s=='quiet', estimated integration errors are not shown  
%
%OUTPUT:
%
%-z: integral. It is a linear vector if x is a matrix.
%
%USAGE:
%
%z=NC_integrate(y,x,n,s)
%
%A. Marinoni 21/02/2012

switch nargin
   case 0
      disp('Insufficient number of input')
      z=[];
      return

   case 1
      disp('Radial grid not provided, assuming dx=1')
      x=1;
      disp('Order not provided, assuming trapezoid')
      n=1;
      y=varargin{1};
      s='show';

   case 2
      disp('Order not provided, assuming trapezoid')
      n=1;
      y=varargin{1};
      x=varargin{2};
      s='show';
      
   case 3
      y=varargin{1};
      x=varargin{2};
      n=varargin{3};
      s='show';
      if or(isempty(n),~isnumeric(n))
         n=1;
	 disp('Invalid order, forcing trapezoid rule')
      end
      if or(n>4,n<1)
         disp('Invalid order, forcing trapezoid rule')
         n=1;
      end
      
   case 4
      y=varargin{1};
      x=varargin{2};
      n=varargin{3};
      s=varargin{4};
      if or(isempty(n),~isnumeric(n))
         n=1;
	 disp('Invalid order, forcing trapezoid rule')
      end
      if or(n>4,n<1)
         disp('Invalid order, forcing trapezoid rule')
         n=1;
      end
end

%Input conditions
if size(y,1)<2
   disp('Nothing to integrate')
   z=[];
   return
end   
if size(x,1)>1
   if and(size(x,2)==1,size(x,1)==size(y,1))
      dx=diff(x,1)*ones(1,size(y,2));
   elseif find(size(x)-size(y))
      disp('X and Y vectors must have the same number of rows and columns')
      z=[];
      return
   else
      dx=diff(x,1);
   end
end
if size(x,1)==1
   if size(x,2)==1
      dx=x;
   elseif size(x,2)==size(y,2)
      dx=ones(size(y,1)-1,1)*x;
   else
      disp('X and Y vectors must have the same number of columns')
      z=[];
      return
   end
end

if ~or(strcmp(s,'quiet'),strcmp(s,'QUIET'))
   if and(n>1,size(x,1)>1)
      disp('WARNING: higher order formulas are implemented only for equally spaced grids.')
      disp('         Use the first order method if grids are not equally spaced')
   end
end

dump=mean(dx,1);
dump2=[(2*dump).^5/2880;(3*dump).^5/6480;(4*dump).^7/1935360]/(dump.^3/12);
dump2=[ones(1,size(y,2));dump2];
if ~or(strcmp(s,'quiet'),strcmp(s,'QUIET'))
   disp(' ')
   disp('Estimated integration errors with respect to first order method:')
   disp(strcat([strvcat(['n=1 ';'n=2 ';'n=3 ';'n=4 ']),strvcat([' ';' ';' ';' ']),num2str(dump2)])) 
end

switch n
   case 1
      z=0.5*sum(dx.*(y(2:end,:)+y(1:end-1,:)));
      
   case 2
      if size(dx,1)==1
         dx=dx*2;
      else
         dx=dx(1:2:end-1,:)+dx(2:2:end,:);
      end
      y=y(1:2:end-2,:)+4*y(2:2:end-1,:)+y(3:2:end,:);
      z=sum(dx.*y,1)/6;
   
   case 3
      disp('WARNING: the 2nd order method is usually better for a given, fixed, grid')
      if size(dx,1)==1
         dx=dx*3;
      else
         dx=dx(1:3:end-2,:)+dx(2:3:end-1,:)+dx(3:3:end,:);
      end
      y=y(1:3:end-3,:)+3*y(2:3:end-2,:)+3*y(3:3:end-1,:)+y(4:3:end,:);
      z=sum(dx.*y,1)/8;
   
   case 4
      if size(dx,1)==1
         dx=dx*4;
      else
         dx=dx(1:4:end-3,:)+dx(2:4:end-2,:)+dx(3:4:end-1,:)+dx(4:4:end,:);
      end
      y=7*y(1:4:end-4,:)+32*y(2:4:end-3,:)+12*y(3:4:end-2,:)+32*y(4:4:end-1,:)+7*y(5:4:end,:);
      z=sum(dx.*y,1)/90;
      
end
