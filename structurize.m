function S = structurize(varargin)
%This function puts all the input arguments into a structure format
%
%INPUT: fields
%
%OUTPUT: name of the desired structure
%
%USAGE: S=strucutrize(input1,input2,...,inputN)
%
%A. Marinoni, 23/08/2012

for k=1:nargin
   eval(strcat(['S.',inputname(k),'=varargin{k};']));
end
