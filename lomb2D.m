function [psd,f,kvec]=lomb2D(yin,varargin)

%This function evaluates the time-space spectrum of a two dimensional
%input signal using the Lomb algorithm extended to handle complex data. 
%The Lomb algorithm is used to handle unevenly sampled data.
%
% INPUTS:
%
%    yin      : input data, it has to be a 2D vector (time,space) [No default]
%    x        : input spatial coordinates [Default, 1,2,3...number of columns of yin]
%    kvec     : spatial angular frequencies to be evaluated [Default, calculated from x vector]
%    nwind    : window to avoid lobe leakage [Default, Hanning with nfft points]
%    nfft     : points over which time-FFT has to be computed [Default, length of window]
%    noverlap : overlapping points between segments [Default, nfft/2]
%    fs       : sampling frequency of data points [Default, 1 Hz]
%
% OUTPUTS:
%
%    psd      : 2D Power Spectral Density (space,time)
%    kvec     : spatial angular frequency [rad/dim(x)]
%    f        : spatial time frequency [rad*dim(fs)]
%
% USAGE: 
%
%    [psd,f,kvec]=lomb(yin,x,k,nwind,nfft,noverlap,fs)
%
%If k is a number, it is understood to be the nember of points from the two Nyquist frequencies 
%If nwind is a number, it is meant as the length of a Hanning window
%
%Data a detrended by subtracting the best-fit straight line from each window
%
%A. Marinoni, 12/07/2011

yout=nan;
fout=nan;
kout=nan;

%Check the input
if nargin<1
   disp('Insufficient number of inputs')
   return
end

[n1,n2]=size(yin);
if or(~isa(yin,'double'),or(n1==1,n2==1))
   disp('The input signal must be double-precision real 2D vector.');
   return
end

defaults={[1:n2],n2,4096,4096,2048,4e6};
variab={'xin','kvec','nwind','nfft','nov','fs'};

lvarar=length(varargin)+1;

for j=1:lvarar-1
   if or(isempty(varargin{j}),~isnumeric(varargin{j}))
      eval(strcat([variab{j},'=[',num2str(defaults{j}),'];']));
   else
      eval(strcat([variab{j},'=[',num2str(varargin{j}),'];']));
   end
end
for j=lvarar:length(variab)
   eval(strcat([variab{j},'=[',num2str(defaults{j})],'];'));
end

kvec
%If kvec is a number, it is assumed to be equal to the number of points
if length(kvec)==1
   dx=mean(diff(xin));
   kvec=pi*linspace(-fix(1/dx),fix(1/dx),kvec);
end
nk=length(kvec);

if length(nwind)==1
   win=hann(nwind);
end

%Calculating the number of columns of the FFT
ncol=fix((n1-nov)/(nwind-nov));

disp(strcat(['Data divided into ',num2str(ncol),' segments of ',num2str((nwind-1)/fs),' s']))

% Pre-process X
xin=xin(:);
colindex=1+(0:(ncol-1))*(nwind-nov);

f=[0:ceil((nfft-1)/2)]/nfft*fs;
f=f(:);

%Hanning window in space
win=win*hann(n2)';
winp=sum(sum(win));
kx=exp(i*xin*kvec);
psd=zeros(ceil(nfft/2),nk);

for j=1:ncol
   
   dump=detrend(yin(colindex(j):colindex(j)+nfft-1,:));
   
   %Windowing
   dump=dump.*win;
   
   %Compute FFT
   dump=fft(dump,nfft);
   dump2=dump(1:ceil(nfft/2),:)*kx;

   %Implementing the modified Lomb algorithm
   psd=psd+dump2.*conj(dump2);

end   

%Normalizations
psd=psd/j/n2/winp;
