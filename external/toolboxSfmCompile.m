function toolboxSfmCompile()
% Compiles all the private routines
%
% assumes located in toolbox root directory
%
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

disp('Compiling.......................................');
savepwd=pwd; cd(fileparts(mfilename('fullpath'))); cd('../');

dir = 'sfm/private/sba/';
dirSBA = [ pwd '/external/sba/'];

disp('************* Compiling SBA');
cd(dirSBA);
delete('sba_levmar.o'); delete('sba_levmar_wrap.o'); delete('sba_lapack.o');
delete('sba_crsm.o'); delete('sba_chkjac.o'); delete('libsba.a');
switch computer
  case 'PCWIN',
    % You need this variable as an environment variable
    %       set INCLUDE="C:\Program Files\Microsoft Visual Studio
    %       8\VC\include"
    system('nmake /f Makefile.vc clean');
    system('nmake /f Makefile.vc sba.lib');
  case {'GLNX86', 'GLNXA64', 'i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
    system('make libsba.a');
end
disp('************* Done compiling SBA');

disp('************* Compiling Matlab SBA');
cd matlab
switch computer
  case 'PCWIN',
    system('nmake /f Makefile.w32 clean');
    system('nmake /f Makefile.w32 sba.mexw32');
  case {'GLNX86', 'GLNXA64', 'i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
    if exist('OCTAVE_VERSION','builtin')
	  mkoctfile --mex ./sba.c -I../ -lsba -L../
    else
	  % Matlab
	  opts={'CXX=g++-4.1' 'CC=g++-4.1' 'LD=g++-4.1' '-I..' '-O' '-l' 'mwlapack' '-l' 'mwblas' '-l' 'dl' '../libsba.a'  '/usr/lib/atlas/libblas.a' '/usr/lib/libgfortran.so.3'};
	  mex( 'sba.c', opts{:} );
	end
end
cd ../../..
disp('************* Done compiling Matlab SBA');

disp('************* Compiling Vincent''s projection routines');
cd sfm/private/sba
switch computer
  case 'PCWIN',
    system('cl /nologo /O2 sbaProjection.c /link /dll /out:sbaProjection.dll');
  case {'GLNX86','GLNXA64','i686-pc-linux-gnu', 'x86_64-pc-linux-gnu'},
    system('gcc -fPIC -O3 -shared -o sbaProjection.so sbaProjection.c');
end
cd ../../..
disp('************* Done compiling Vincent''s projection routines');

disp('************* Compiling computeCsfmInfimum, refineExteriorOrientation and msfmRotationDerivative');

rd=fileparts(mfilename('fullpath')); rd=rd(1:end-9);

% general compile options (can make architecture specific)
if exist('OCTAVE_VERSION','builtin') opts = {'-o'};
else opts = {'-output'};
  % if you get warnings on linux, you couldforce the gcc version this way
  % opts = {'CXX=g++-4.1' 'CC=g++-4.1' 'LD=g++-4.1' '-l' ...
  % 'mwlapack' '-l' 'mwblas' '-output' };
end

% general compile options
fs={'computeCsfmInfimumMex','refineExteriorOrientationMex', ...
  'msfmRotationDerivativeMex'};
ds={'nrsfm', 'sfm', 'msfm'};
for i=1:length(fs), mex([rd '/' ds{i} '/private/' fs{i} '.c'],...
    opts{:},[rd '/' ds{i} '/private/' fs{i} '.' mexext]); end

cd(savepwd); disp('..................................Done Compiling');

end
