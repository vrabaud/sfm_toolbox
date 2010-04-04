function toolboxSfmCompile()
% Compiles all the private routines
%
% assumes located in toolbox root directory
%
% if you get warnings on linux because yourgcc is too old/new, you can
% force the gcc version this way
%  opts = {'CXX=g++-4.1' 'CC=g++-4.1' 'LD=g++-4.1' etc
% USAGE
%  toolboxCompile
%
% INPUTS
%
% OUTPUTS
%
% EXAMPLE
%
% See also
%
% Done

if strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64')
  gcc_ver = '4.2';
  gcc_extra = [ 'CXX=g++-' gcc_ver ' CC=g++-' gcc_ver ' LD=g++-' gcc_ver ];
end

disp('Compiling.......................................');
savepwd=pwd; cd(fileparts(mfilename('fullpath'))); cd('../');

dirSBA = [ pwd '/external/sba/'];

cd(dirSBA);

% delete previous object files
fileList = {'sba_levmar.o', 'sba_levmar_wrap.o', 'sba_lapack.o', ...
  'sba_crsm.o', 'sba_chkjac.o', 'libsba.a', 'sba.lib', ...
  'matlab/sba.mexw32'};
for i=1:length(fileList)
  if exist(fileList{i}, 'file'), delete(fileList{i}); end
end

switch computer
  case 'PCWIN',
    % Matlab on Windows 32 (not sure about 64)
    % You need this variable as an environment variable
    %       set INCLUDE="C:\Program Files\Microsoft Visual Studio
    %       8\VC\include"
    system('nmake /f Makefile.vc sba.lib');
  case {'GLNX86', 'GLNXA64', 'i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
    % Matlab and Octave on Linux
    system('make libsba.a');
end

cd matlab
switch computer
  case 'PCWIN',
    % Matlab on Windows 32 (not sure about 64)
    system('nmake /f Makefile.w32 sba.mexw32');
  case {'GLNX86'},
	% Matlab on Linux
	eval([ 'mex ' gcc_extra '-I.. -O -ldl sba.c ../libsba.a ' ...
	  '/usr/lib/liblapack.a /usr/lib/atlas/libblas.a ' ...
	  '/usr/lib/libgfortran.so.3']);
  case {'GLNXA64'},
	% Matlab on Linux 64
	eval([ 'mex -lgfortran ' gcc_extra ' -I.. -O -ldl  sba.c ' ...
	  '../libsba.a /usr/lib/liblapack_pic.a ' ...
	  '../../blasLapack/linux64/libf77blas.a ' ...
	  '../../blasLapack/linux64/libatlas.a '
	  ]);
  case {'i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
	% Octave on Linux
    mkoctfile --mex ./sba.c -I../ -lsba -L../
end

cd ../../..

cd sfm/private/sba
switch computer
  case 'PCWIN',
    % Matlab on Windows 32 (not sure about 64)
    system('cl /nologo /O2 sbaProjection.c /link /dll /out:sbaProjection.dll');
  case {'GLNX86','GLNXA64','i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
    % Matlab and Octave on Linux
    system('gcc -Wall -fPIC -O3 -shared -o sbaProjection.so sbaProjection.c');
end
cd ../../..

rd=fileparts(mfilename('fullpath')); rd=rd(1:end-9);

% general compile options (can make architecture specific)
optsAfter={};
switch computer
  case 'PCWIN',
    % Matlab on Windows 32 (not sure about 64)
	lapacklib = fullfile(matlabroot, 'extern', 'lib', 'win32', ...
	  'microsoft', 'libmwlapack.lib');
	blaslib = fullfile(matlabroot, ...
	  'extern', 'lib', 'win32', 'microsoft', 'libmwblas.lib');
	opts={'-output'};
	optsAfter = {lapacklib, blaslib};
  case {'GLNX86', 'GLNXA64'},
    % Matlab on Linux
	% if you get warnings on linux, you could force the gcc version by
	% adding those options: 'CXX=g++-4.1' 'CC=g++-4.1' 'LD=g++-4.1'
	opts = { ['CXX=g++-' gcc_ver] ['CC=g++-' gcc_ver] ['LD=g++-' ...
	  gcc_ver] '-largeArrayDims' '-l' 'mwlapack' '-l' 'mwblas' '-output'};
  case {'i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
	% Octave on Linux
	opts = {'-o'};
end

% general compile options
fs={'computeCsfmInfimumMex','refineExteriorOrientationMex', ...
  'msfmRotationDerivativeMex'};
ds={'nrsfm', 'sfm', 'msfm'};
for i=1:length(fs)
  mex(opts{:},[rd '/' ds{i} '/private/' fs{i} '.' mexext], ...
    [rd '/' ds{i} '/private/' fs{i} '.c'], optsAfter{:});
end

cd(savepwd); disp('..................................Done Compiling');

end
