function gen_onetwo_netcdf(filename)

%This function unzips the TRPLTFIL onetwo
%output file given in input and calls
%the fortran routine prplt to convert into 
%a netcdf format file with the same name as filename.
%
%INPUT:
%
%- filename: name of the gunzipped trpltfil file
%
%OUTPUT:
%
%- the netcdf file filename.nc is saved in the same directory.
%
%USAGE:
% 
%>> gen_onetwo_netcdf(filename);
%
%A. Marinoni, 20/06/2012

if nargin<1
   filename=input('Provide input gunzipped filename as a string');
end
if ~strcmp(filename(end-2:end),'.gz')
   filename=strcat(filename,'.gz');
   disp('Added .gz extension to filename')
end
if exist(filename)~=2
   disp('The input file does not exist')
   return
end

disp(strcat(['Unzipping ',filename,' ...']))
[s,w]=unix(strcat(['gzip -d ',filename]));
if s
   disp(strcat(['Problem with unix command, ',w]))
end
disp('Generating input for preplt routine ...')
[s,w]=unix(strcat(['cp ',filename(1:end-3),' trpltfil']));
if s
   disp(strcat(['Problem with unix command cp: ',w]))
end
disp('Running preplt routine ...')
[s,w]=unix('preplt');
if s
   disp(strcat(['Problem with unix command preplt: ',w]))
end
disp('Renaming output files ...')
[s,w]=unix(strcat(['mv trpltout.nc ',filename(1:5),'out',filename(9:end-3),'.nc']));
if s
   disp(strcat(['Problem with unix command mv: ',w]))
else
   disp(strcat(['Netcdf file ',filename(1:5),'out',filename(9:end-3),'.nc was generated']))
end
[s,w]=unix(strcat(['rm trpltfil trpltout']));
if s
   disp(strcat(['Problem with unix command rm: ',w]))
end
disp(strcat(['Gunzipping ',filename]))
[s,w]=unix(strcat(['gzip ',filename(1:end-3)]));
if s
   disp(strcat(['Problem with unix command gzip: ',w]))
else
   disp(strcat([filename(1:end-3),' has been gzipped again']))
end
