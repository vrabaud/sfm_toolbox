function scriptName = gcLaunch( nitr, matlabCommand, exePath, savePath,...
  scriptPath )

temp = clock; scriptName = [ 'script-' int2str(temp(4)) '-' ...
  int2str(temp(5)) '-' int2str(temp(6)) ];

% Create the script file
fid=fopen('launchJob.sh','w');
fprintf(fid,'#!/bin/bash \n');
fprintf(fid,'nitr=%d\n', nitr );
fprintf(fid,'matlabCommand=''%s''\n', matlabCommand);
fprintf(fid,'executePath=''%s''\n',exePath);
fprintf(fid,'matlabPath=''/projects/p1/software/matlab2007b/bin/matlab -nodesktop''\n');
fprintf(fid,'\n');
%fprintf(fid,'# kill old jobs and remove temp variables\n');
%fprintf(fid,'qdel -u vrabaud\n');
fprintf(fid,'find %s -name ''script*'' | xargs rm\n', savePath );
fprintf(fid,'find %s -name ''script-*.o*'' | xargs rm\n', exePath );
fprintf(fid,'find %s -name ''script-*.e*'' | xargs rm\n', exePath );
%fprintf(fid,'find ~/temp -name ''temp*'' | xargs rm\n');
fprintf(fid,'\n');
fprintf(fid,'# For each script\n');
fprintf(fid,'i=1\n');
fprintf(fid,['for clus in  '...
  'renderx01 renderx02 renderx03 renderx04 renderx05 renderx06 renderx07' ...
  'renderx08 renderx09 renderx10 renderx12 renderx13 renderx14 ' ...
  'renderx15 rendera01 rendera02 rendera03 rendera04 rendera05 ' ...
  'rendera06 rendera07 rendera08 rendera10 rendera11 rendera12 ' ...
  'rendera13 rendera14 rendera16\n']);
fprintf(fid,'do\n');
fprintf(fid,'  # Create the script\n');
fprintf(fid,['  todo=`ssh -q -n vrabaud@${clus} "cd ${executePath}; ${matlabPath} -r \\"$' ...
  '{matlabCommand}(${i},''%s'',''%s'')\\" " ` &\n'], scriptName,savePath);
fprintf(fid,'  \n');
fprintf(fid,'  # If no possible allocation, wait and resubmit\n');
fprintf(fid,'  while true\n');
fprintf(fid,'  do\n');
fprintf(fid,'    # Check the number of running MATLAB licenses\n');
fprintf(fid,'    nLicense=$(/projects/p1/software/matlab2007b/etc/lmstat -f MATLAB | grep -o "Total of [0-9]* licenses in use" | grep -o "[0-9][0-9]*")\n');
fprintf(fid,'    \n');
fprintf(fid,'    if [ $nLicense -gt 490 ]\n');
fprintf(fid,'    then\n');
fprintf(fid,'      sleep 1\n');
fprintf(fid,'      continue\n');
fprintf(fid,'    fi\n');
fprintf(fid,'    \n');
fprintf(fid,'    $todo\n');
fprintf(fid,'    break\n');
fprintf(fid,'  done\n');
fprintf(fid,'  ((i=i+1))\n');
fprintf(fid,'done\n');

fclose(fid);

% Kill any running job on the cluster
%system('ssh -q -n vrabaud@fwg-cs0 "kill -9 -1"');
% Copy the script on the cluster
system(['mv launchJob.sh ' exePath]);
system(['chmod 755 ' exePath '/launchJob.sh']);
system([exePath '/./launchJob.sh &']);


