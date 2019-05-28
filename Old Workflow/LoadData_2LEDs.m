function [m,data] = LoadData_2LEDs(m)
tic
m.nLEDs = 2;
runnum=1;
runtoload=[m.CCDpath{m.ttu(runnum)} '/' m.stimlist{m.ttu(runnum)} '/'];
load(cell2mat(fullfile(m.pathofmetadata,m.metadataruns(m.ttu(runnum)))))

m.fr = frameRate; m.binning = binSize;
m.nFrames=length(dir(runtoload))-2;
m.t = linspace(0,(tpre+tpost+tstim),floor(nFrames/2));
m.height = height; m.width = width;
m.tpre = tpre; m.tstim = tstim; m.tpost = tpost;

delete(gcp('nocreate'));
all_data=zeros(height,width,nFrames);
blue_avg=zeros(height,width,floor(nFrames/2));
green_avg=zeros(height,width,floor(nFrames/2));
summer=0;
h = waitbar(0,'Averaging...');
for runnum=1:length(m.ttu)
    bla=[m.CCDpath{m.ttu(runnum)} '/' m.stimlist{m.ttu(runnum)} '/' m.stimlist{m.ttu(runnum)}(1:4)];
    for j=0:m.nFrames-1
        fid = fopen([bla repmat('0',[1,10-size(num2str(j),2)]) num2str(j) '.dat'],'r','l');
        all_data(:,:,j+1) = fread(fid,[m.height m.width],'uint16','l');
        fclose(fid);
    end
    blue_avg=blue_avg+all_data(:,:,1:2:end-1);
    green_avg=green_avg+all_data(:,:,2:2:end);
    summer=summer+1;
    waitbar(runnum/numel(m.ttu),h,['finished ' mat2str(runnum) ' out of ' mat2str(length(m.ttu)) ' runs.']);
end
close(h)
data.b=rot90(blue_avg/summer,m.nrot);
data.g=rot90(green_avg/summer,m.nrot);
disp(['Total time - ' mat2str(toc) ' seconds, taking about ' mat2str(toc/numel(m.ttu)) ' seconds per trial'])

