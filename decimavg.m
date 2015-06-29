function [outdata,outtime]=decimavg(indata,intime,ofreq,varargin)

%This function downsamples the input data to a new sampling 
%frequency chosen by the user; the algorithm will then decimate
%the signal to the closest sampling frequency. 
%Before downsampling, it averages the input signal over points 
%corresponding to higher frequencies than the resulting Nyquist 
%frequency of the output.
%The average is perfomed with a Hann window.
%The filter used in the downsampling process is, by default, a type I
%Chebycheff filter of 9th order. The user can select other filters as specified below.
%For matrices, it filters along columns
%
%The function assumes 1D vectors
%
%USAGE:
% 
%[outdata,outtime]=decimavg(indata,intime,ofreq,filter,filter order)
%
%outdata = downsampled data
%indata = original data
%intime = original time dataset, the routine detects the sampling frequency 
%If the array has unitary length, then the number provided is assumed to be the sampling frequency 
%and it is assumed to start at t=0
%outtime = resulting time start. 
%ofreq = resulting sampling rate in Hz  
%
%OPTIONAL INPUT FILTER: the filter can be chosen among the following
%'cheby1': Chebycheff type I filter
%'cheby2': Chebycheff type II filter
%'elliptical': Elliptical filter
%'buttord': Butterworth filter
%
%OPTIONAL INPUT FILTER ORDER: it can be anything larger than zero
%
%A. Marinoni 30/06/2011  

%Checks on the input
if nargin<3
   disp('Insufficient number of inputs')
   outdata=nan;
   outtime=nan;
   return
end

nt=length(intime);
[n1,n2]=size(indata);
if and(nt>1,nt~=n1)
   disp('Not matching data and time: input data have different lengths')
   outdata=nan;
   outtime=nan;
   return
end

if isempty(indata) || isempty(intime) || ~isa(indata,'double') || ~isa(intime,'double')
    disp('The input signal must be double-precision vectors.');
    outdata=nan;
    outtime=nan;
    return
end

%Sampling rate of inputs
if nt==1
   ifreq=intime;
else
   ifreq=1/mean(diff(intime));
   disp(strcat(['Detected frequency: ',num2str(ifreq),' [in inverse units of input time]']))
end

if (ofreq > 0.5*ifreq ) || (ofreq <= 0)
    disp('The chosen new sampling frequency must be positive and less than half of the original sampling frequency');
    outdata=nan;
    outtime=nan;
    return
end

if ofreq==ifreq
   disp('Unitary downsampling factor. Output = input')
   outdata=indata;
   outtime=intime;
   return
end

if length(varargin)<2
   varargin{2}=[];
end

filtype=varargin{1};
if isempty(filtype)
   filtype='cheby1';
end

nfilt=round(varargin{2});
if isempty(nfilt)
   nfilt=9;
end

%Column vectors
%indata=indata(:);
intime=intime(:);

%Averaging to reduce noise

%windows length
wl=floor(ifreq/ofreq);

%window choice, after a few tests, the Hann function is  
%set as default

wind=hann(wl);
%wind=ones(wl,1);
%wind=hamming(wl);
%wind=gausswin(wl);
%wind=chebwin(wl);

%Choosing the filter 

switch filtype

   case 'cheby1'

      rip = 0.05;	% passband ripple in dB
      [b,a] = cheby1(nfilt, rip, 1/wl);
      while all(b==0) || (abs(filtmag_db(b,a,1/wl)+rip)>1e-6)
         nfilt = nfilt - 1;
       	 disp('Null numerator or too much ripple at cut-off, decreasing filter order')
	 if nfilt == 0
            break
         end
         [b,a] = cheby1(nfilt, rip, 1/wl);
      end
      disp(strcat(['Using ',num2str(nfilt),'th Chebychef type I filter']))
      
   case 'buttord'
      
      [b,a] = butter(nfilt, 1/wl);
      while all(b==0) 
         nfilt = nfilt - 1;
       	 disp('Null numerator, decreasing filter order')
         if nfilt == 0
            break
         end
         [b,a] = butter(nfilt, 1/wl);
      end
      disp(strcat(['Using ',num2str(nfilt),'th Butterworth filter']))
      
   case 'cheby2'
   
      rip = 20;	% stopband ripple in dB
      [b,a] = cheby2(nfilt, rip, 1/wl);
      while all(b==0) || (abs(filtmag_db(b,a,1/wl)+rip)>1e-6)
         nfilt = nfilt - 1;
	 disp('Null numerator. Decreasing filter order')
         if nfilt == 0
            break
         end
         [b,a] = cheby2(nfilt, rip, 1/wl);
      end
      disp(strcat(['Using ',num2str(nfilt),'th Chebychef type II filter']))
      
   case 'elliptical'
   
      rip = 0.05;	% passband ripple in dB
      att = 20;        % stopband attenuation in dB
      [b,a] = ellip(nfilt, rip, att, 1/wl);
      while all(b==0) || (abs(filtmag_db(b,a,1/wl)+rip)>1e-6)
         nfilt = nfilt - 1;
       	 disp('Null numerator or too much ripple at cut-off, decreasing filter order')
         if nfilt == 0
            break
         end
         [b,a] = ellip(nfilt, rip, att, 1/wl);
      end
      disp(strcat(['Using ',num2str(nfilt),'th Elliptic filter']))
      
   otherwise
   
      disp('Unrecognized filter type. Returning')
      outdata=indata;
      outtime=intime;
      return

end
if nfilt == 0
   error('Bad filter design, likely R is too big; try mult. decimation (R=R1*R2).')
end

%Convolving before filtering
[n1,n2]=size(indata);
outdata=[];

for i=1:n2

   sig=conv(wind,indata(:,i));
   sig=sig/sum(wind);
   %To avoid partial sums at end points in the convolution
   sig=sig(wl:end-wl+1);
   M=length(sig);

   %Filtering, (making use of on MATLAB function "decimate"
   %after having added the possibility of using a other types of filters

   % be sure to filter in both directions to make sure the filtered data has zero phase
   % make a data vector properly pre- and ap- pended to filter forwards and back
   % so end effects can be obliterated.
   dump = filtfilt(b,a,sig);
   outdata = [outdata dump(1:wl:M)];
   
end    
%Generating time data

if nt>1
   outtime=intime(ceil(wl/2):wl:n1-floor(wl/2));
else
   outtime=linspace((ceil(wl/2)-1)/ifreq,(n1-floor(wl/2)-1)/ifreq,size(outdata,1));
   outtime=outtime(:);
end
%--------------------------------------------------------------------------
function H = filtmag_db(b,a,f)
%FILTMAG_DB Find filter's magnitude response in decibels at given frequency.

nb = length(b);
na = length(a);
top = exp(-j*(0:nb-1)*pi*f)*b(:);
bot = exp(-j*(0:na-1)*pi*f)*a(:);

H = 20*log10(abs(top/bot));


%--------------------------------------------------------------------------

%Observations on filters:
%
%A spectrum consisting of one hundred harmonics (randomly chosen in amplitude 
%and frequency) was generated along with white gaussian noise at 15 dB S/R.
%Resulting FFT spectra were compared. Digitizing frequency was 10 MHz and
%downsampled to 1 MHz.
%The minimum order of the Elliptical filter was 7th, all filters were
%therefore compared at this order.
%
%The Butterworth filter tends to perform worse in recovering peaks at 
%frequencies close to the Nyquist frequency of the resampled dataset.
%The loss of gain acts on the noise level by lowering it compared to other
%filters. The loss was quantified as about a factor 3.
%
%Chebychef type II filter starts cutting the signal at about 440 kHz.
%
%Chebycheff type I and Elliptical fiters were found to be approximately 
%equivalent along the whole frequency spectrum. Since the Chebycheff type I
%filter can be built with a larger order for the same decimation factor,
%this was chosen as default filter.
%
%
%Obeservations on averaging:
%
%Several windows were tested: Rectangular, Gaussian, Chebycheff, Hann and Hamming.
%The rectangular window cuts the signal and the noise more that all the others.
%The difference being larger and larger at higher frequencies.
%The other windows are pretty much equivalent. 
%The S/N ratio is approximately constant for all the windows.
%The Hann window was chosen as default. 
