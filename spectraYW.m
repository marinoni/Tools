function [P,varargout] = spectraYW(varargin)

%This function computes the 2D spectrum
%of an input signal using standard FFT
%for time (along dimension 1), and the 
%Maximum Entropy Method implemented with the 
%Yule-Walker equations (along dimension 2).
%Data are linear detrended and time windowed with a 
%Hanning(NFFT) window
%
%INPUT:
%
% X         : signal 
% NFFT      : number of points to compute the FFT in time.
%             Default: 2048
% OV        : Number of overlapping points between windows in temporal
%             FFT. Default: 1024
% N         : order of the Autoregressive model used to 
%             estimate the spatial spectrum. Default: 4
% tsampfreq : time sampling. Default: 1 s
% ssampfreq : spatial sampling. Default: 1 cm
%
%OUTPUT:
%
% P : Spectrum (complex in general)
% f : time frequency in units of [1/tsampfreq]. To be multiplied by 
%     2*pi in order to obtain units of rad/s
% k : spatial frequency (to be multiplied by 2*pi in order
%     to obtain rad/cm
%
%USAGE:
%
% [P,f,k]=spectraYW(x,NFFT,OV,N,tsampfreq,ssampfreq);
%
% A. Marinoni, 29/03/2012

def={2048;1024;4;1;1};
variab={'x';'nfft';'ov';'N';'tsampfreq';'ssampfreq'};

switch nargin

   case 0
      disp('Insufficient number of inputs')
      P=[];
      return

   otherwise
      x=varargin{1};
      for i=2:nargin
         if isempty(varargin{i})
	    eval(strcat([variab{i},'=[];']));
	 else
            eval(strcat([variab{i},'=',num2str(varargin{i}),';']));
         end
      end
      for i=nargin+1:length(variab)-1
         if isempty(def{i-1})
	    eval(strcat([variab{i},'=[];']));
	 else
	    eval(strcat([variab{i},'=',num2str(def{i-1})],';'));
	 end
      end
      
end

%Checking inputs
for i=1:length(variab)
   if ~eval(strcat(['isnumeric(',variab{i},')']))
      disp('Incorrect format for ',variab{i})
      P=[];
      return
   end
end
if ndims(x)~=2
   disp('Input signal must be two dimensional')
   P=[];
   return;
end
if ~isreal(x)
   disp('Input signal must be real')
   P=[];
   return;
end
if ov>nfft
   disp('Overlap cannot exceed time window length. Forcing nov=nfft/2')
   ov=floor(nfft/2);
end      
[n1,n2]=size(x);
%Calculating the number of columns of the FFT
ncol=fix((n1-ov)/(nfft-ov));
disp(strcat(['Data divided into ',num2str(ncol),...
             ' segments of ',num2str((nfft-1)/tsampfreq),' s']))

% Pre-process X
rowindex=1+(0:(ncol-1))*(nfft-ov);
colindex=(1:nfft)';

x=x.*(ones(n1,1)*hann(n2)');
y=[];
disp('Computing temporal FFT')
for j=1:ncol
   y=[y x(rowindex(j):rowindex(j)+nfft-1,:)]; 
end

%Detrend
%y=detrend(y);

%Windowing
win=hann(nfft);
y=y.*(win*ones(1,ncol*n2));
  
%Compute FFT
z=fft(y,nfft);
z=z(1:ceil(nfft/2),:);

y=0;
for i=1:ncol
   y=z(:,1+(j-1)*n2:j*n2)+y;
end
y=y/ncol;

%Time frequency vector
f=[0:nfft-1];
f=f/(nfft-1)*tsampfreq;
f=f(1:ceil(nfft/2));
k=[0:n2-1]/(n2-1)*ssampfreq;
k=k-k(ceil(n2/2));

nf=length(f);

disp('Computing spatial Fourier transform')
%Computing Spatial spectrum and averaging over windows
P=zeros(n2,nf);
for i=1:nf
   P(:,i)=pyulear(y(i,:),N,n2,'twosided');
end

%Transposing and averaging
P=P'/ncol;

%Renormalization to compensate for the power of the windows
fac=sum(win)/nfft;
fac2=sum(hann(n2))/n2;
P=P/fac^2/fac2^2;
P=fftshift(P,2);
varargout{1}=f;
varargout{2}=k;
