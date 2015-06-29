function [y,varargout]=bspline_rec(varargin)

%This function reconstructs a function fitted with B-splines
%of order n, coefficents c and knot locations kl, on a vector xs
%
%INPUT:
%
%-  c: spline coefficients
%- kl: knot locations [Default: linear array from 0 to 1 of length(c) elements]
%-  n: spline order [Default: cubic]
%-  x: linear array on which the function is to be evaluated [Default: linear 101 elements in [0,1] ]
%      If xs is a one digit integer, a linear array in [0,1] with length=xs is assumed
%
%OUTPUT:
%
%- y: interpolated function
%- x: array on which y is evaluated [Optional]
%
%USAGE:
%
% [y,x]=bspline_rec(c,kl,n,xs)
%
%A. Marinoni 06/03/2012


switch nargin
   case 0
      disp('Insufficient number of input')
      y=[];
      x=[];
      return

   case 1
      c=varargin{1};
      if or(isempty(c),~isnumeric(c))
          disp('The format of coefficients is not correct')
	  varargout{1}=[];
	  y=[];
	  return
      end
      n=3;
      disp(strcat(['Knot locations not correct: assuming linear vector [0,1] of length ',num2str(length(c)+2*(n-1))]))
      kl=linspace(0,1,length(c));
      disp('Spline order not provided, assuming cubic')
      disp('Radial grid not provided, assuming 101 elements linear vector [0,1]')
      x=linspace(0,1,101);
      %Add n+1 additional knots to provide full basis of independent splines 
      kl=[zeros(1,ceil((n+1)/2)) kl ones(1,floor((n+1)/2))];

   case 2
      c=varargin{1};
      kl=varargin{2};
      disp('Spline order not provided, assuming cubic')
      n=3;
      disp('Radial grid not provided, assuming 101 elements linear vector [0,1]')
      x=linspace(0,1,101);
      
   case 3
      c=varargin{1};
      kl=varargin{2};
      n=varargin{3};
      disp('Radial grid not provided, assuming 101 elements linear vector [0,1]')
      x=linspace(0,1,101);
      
   otherwise
      c=varargin{1};
      kl=varargin{2};
      n=varargin{3};
      x=varargin{4};
      
end

%Conditions on inputs
if isstr(c)
   disp('Knots coefficients cannot be string-type. Aborting')
   y=nan;
   varargout{1}=nan;
   return;
end
if or(isempty(n),~isnumeric(n))
   disp('Incorrect spline order: assuming cubic');
   n=3;
else
   n=int16(n);
end
if n<1
   disp('The order must be >= 1. Forcing cubic order')
   n=3;
end

if or(isempty(kl),~isnumeric(kl))
   disp(strcat(['Knot locations not correct: assuming linear vector [0,1] of length ',num2str(length(c)+2*(n-1))]))
   kl=linspace(0,1,length(c));
   %Add n+1 additional knots to provide full basis of independent splines 
   kl=[zeros(1,ceil((n+1)/2)) kl ones(1,floor((n+1)/2))];
end

if or(isempty(x),~isnumeric(x))
   disp('Incorrect radial grid: assuming 101 elements linear vector [0,1]');
   x=linspace(0,1,101);
end

if length(x)==1
   disp(strcat(['Assuming ',num2str(x),' elements linear vector [0,1]']))
   x=linspace(0,1,int8(x));
end

if length(c)-(length(kl)-(n+1))
   disp('The number of coefficients and knot locations are inconsistent with each other')
   varargout{1}=nan;
   y=nan;
   return
end
m=length(kl);
xl=length(x);
%Making coefficients a row vector
c=c(:)';
%Making x vector a column vector
x=x(:);

%Computing basis functions through Cox-de Boor algorithm

%Initializing
B(1).B=zeros(xl,m-1);
for j=1:m-1
   indx1=find(x>=kl(j),1,'first');
   indx2=find(x<kl(j+1),1,'last');
   B(1).B(indx1:indx2,j)=1;
end 

%Iterating to find higher order basis functions
for h=2:n+1
   B(h).B=zeros(xl,m-1);
   for j=1:m-n-1
      a=kl(j+h-1:j+h)-kl(j:j+1);
      dump=find(a);
      b=zeros(xl,2);
      for i=1:length(dump)
         num=[x-kl(j) kl(j+h)-x];
         b(:,dump(i))=num(:,dump(i))./a(dump(i));
      end
      B(h).B(:,j)=B(h-1).B(:,j).*b(:,1)+B(h-1).B(:,j+1).*b(:,2);
   end
end

%Reconstructing the function
c=ones(xl,1)*c;
y=sum(B(end).B(:,1:m-n-1).*c,2);
varargout{1}=x;
