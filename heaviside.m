function out=heaviside(val)

%Heaviside function. For complex values it returns
%the step function of the real part.
%In zero the step function is defined to return 0.5
%
%A. Marinoni, 16/01/2013

if ~isnumeric(val)
   disp('Input must be a numeric value')
   return
end

val=real(val);
ind=find(val>0);
out(ind)=1;
ind=find(val<0);
out(ind)=0;
ind=find(val==0);
out(ind)=0.5;
