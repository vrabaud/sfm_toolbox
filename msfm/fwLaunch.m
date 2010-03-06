function scriptName = fwLaunch( nitr, matlabCommand, exePath, savePath,...
  scriptPath )

temp = clock; scriptName = [ 'script-' int2str(temp(4)) '-' ...
  int2str(temp(5)) '-' int2str(temp(6)) ];

% Create the script file
fid=fopen('launchJob.sh','w');
fprintf(fid,'#!/bin/bash \n');
fprintf(fid,'nitr=%d\n', nitr );
fprintf(fid,'matlabCommand=''%s''\n', matlabCommand);
fprintf(fid,'executePath=''%s''\n',exePath);
fprintf(fid,'matlabPath=''/apps/matlab/bin/matlab -nosplash -nodesktop''\n');
fprintf(fid,'\n');
%fprintf(fid,'# kill old jobs and remove temp variables\n');
%fprintf(fid,'qdel -u vrabaud\n');
fprintf(fid,'find %s -name ''*'' | xargs rm\n', savePath );
fprintf(fid,'find %s -name ''script-*.o*'' | xargs rm\n', exePath );
fprintf(fid,'find %s -name ''script-*.e*'' | xargs rm\n', exePath );
%fprintf(fid,'find ~/temp -name ''temp*'' | xargs rm\n');
fprintf(fid,'\n');
fprintf(fid,'# For each script\n');
fprintf(fid,'i=1\n');
fprintf(fid,'while [ $i -le $nitr ]\n');
fprintf(fid,'do\n');
fprintf(fid,'  # Create the script\n');
fprintf(fid,'  name=%s/%s-$i\n', scriptPath, scriptName);
fprintf(fid,['  echo "cd ${executePath}; ${matlabPath} -r \\"$' ...
  '{matlabCommand}(${i},''%s'',''%s'')\\" " >> $name\n'], scriptName,savePath);
fprintf(fid,'  chmod 755 $name\n');
fprintf(fid,'  \n');
fprintf(fid,'  # If no possible allocation, wait and resubmit\n');
fprintf(fid,'  while true\n');
fprintf(fid,'  do\n');
fprintf(fid,'    # Check the number of running MATLAB licenses\n');
fprintf(fid,'    nLicense=$(/apps/matlab/etc/lmstat -f MATLAB | grep -o "Total of [0-9]* licenses in use" | grep -o "[0-9][0-9]*")\n');
fprintf(fid,'    \n');
fprintf(fid,'    if [ $nLicense -gt 490 ]\n');
fprintf(fid,'    then\n');
fprintf(fid,'      sleep 1\n');
fprintf(fid,'      continue\n');
fprintf(fid,'    fi\n');
fprintf(fid,'    \n');
fprintf(fid,'    # Check if the script can be submitted properly\n');
fprintf(fid,'    test=`qsub -cwd $name 2>&1 | grep -c submitted`\n');
fprintf(fid,'    if [ $test -gt 0 ]\n');
fprintf(fid,'    then\n');
fprintf(fid,'      break\n');
fprintf(fid,'    fi\n');
fprintf(fid,'    \n');
fprintf(fid,'    sleep 1\n');
fprintf(fid,'  done\n');
fprintf(fid,'  ((i=i+1))\n');
fprintf(fid,'done\n');

fclose(fid);

% Kill any running job on the cluster
%system('ssh -q -n vrabaud@fwg-cs0 "kill -9 -1"');
% Copy the script on the cluster
system(['scp launchJob.sh vrabaud@137.110.131.21:' exePath]);
system( 'rm launchJob.sh' );
% Make it executable
system(['ssh -q -n vrabaud@fwg-cs0 "chmod 755 ' exePath '/launchJob.sh"']);
% Touch it 
system(['ssh -q -n vrabaud@fwg-cs0 "touch ' exePath '/launchJob.sh"']);
% Launch it
system(['ssh -q -n vrabaud@fwg-cs0 "cd ' exePath '; ./launchJob.sh" &']);

% Reconnect to the ssh shares
system('sshfs vrabaud@137.110.131.21:/home/vrabaud /projects/fwHome');
system(['sshfs vrabaud@137.110.131.21:' ...
  '/scratch/vrabaud/ /projects/fwScratch']);
