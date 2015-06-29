function y=gramschmidt(x)

%This function implements the modified Gram-Schimdt decomposition of 
%input matrix x, 
%The routine operates column-wise
%
%INPUT:
%
% - x : input matrix
%
%OUTPUT
%
% - y : output orthogonalized matrix
%
%USAGE:
%
% y=gramschmidt(x)
%
%A. Marinoni, 14/03/2013

y=nan;
if isempty(x)
   disp('Error: empty input');
   return
end
if ~isnumeric(x)
   disp('Error: input values are not numeric')
   return
end
[n,m]=size(x);
y=[];
for i=1:m
   v=x(:,i);
   h=v;
   for j=1:i-1
      h=h-dot(h,y(:,j))/norm(y(:,j))^2*y(:,j);   
   end
   y(:,i)=h;
end
