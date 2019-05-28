function [m,data,datads] = LoadData(m)
tic
runnum=1;
runtoload=[m.CCDpath{m.ttu(runnum)} '/' m.stimlist{m.ttu(runnum)} '/'];
load(cell2mat(fullfile(m.pathofmetadata,m.metadataruns(m.ttu(runnum)))))

m.fr = frameRate; m.binning = binSize;
m.nFrames=length(dir(runtoload))-2;
m.t = linspace(0,(tpre+tpost+tstim),floor(nFrames/3));
m.height = height; m.width = width;
m.tpre = tpre; m.tstim = tstim; m.tpost = tpost;

delete(gcp('nocreate'));
all_data=zeros(height,width,nFrames);

blue_avg=zeros(height,width,floor(nFrames/3));
green_avg=zeros(height,width,floor(nFrames/3));
red_avg=zeros(height,width,floor(nFrames/3));
summer=0;
h = waitbar(0,'Averaging...');
for runnum=1:length(m.ttu)
    bla=[m.CCDpath{m.ttu(runnum)} '/' m.stimlist{m.ttu(runnum)} '/' m.stimlist{m.ttu(runnum)}(1:4)];
    for j=0:m.nFrames-1
        fid = fopen([bla repmat('0',[1,10-size(num2str(j),2)]) num2str(j) '.dat'],'r','l');
        all_data(:,:,j+1) = fread(fid,[m.height m.width],'uint16','l');
        fclose(fid);
    end
    blue_avg=blue_avg+all_data(:,:,1:3:end-1);
    green_avg=green_avg+all_data(:,:,2:3:end);
    red_avg=red_avg+all_data(:,:,3:3:end);
    datads.b(:,:,:,runnum) = imresize(all_data(:,:,1:3:end-1),[32 32]);
    datads.g(:,:,:,runnum) = imresize(all_data(:,:,2:3:end-1),[32 32]);
    datads.r(:,:,:,runnum) = imresize(all_data(:,:,3:3:end-1),[32 32]);

    % If you want TCs for all runs, uncomment this code!
%     if runnum == 1
%         [chbo_tmp,chbr_tmp,chbt_tmp] = convert_mariel(green_avg,red_avg,'g','r',m.interval,534);
%         gcamp_tmp=GcampMcCorrection(blue_avg,chbr_tmp,chbo_tmp,m.interval,m.mua1,m.mua2);
%         m = dda(chbo_tmp,m,'chbo');
%         f = warndlg('Press ok when done picking ROI','Cool?');
%         waitfor(f); 
%         
%     end
%     datac.b(:,:,:,runnum) = all_data(m.y(1):m.y(2),m.x(1):m.x(2),1:3:end-1);
%     datac.g(:,:,:,runnum) = all_data(m.y(1):m.y(2),m.x(1):m.x(2),2:3:end-1);
%     datac.r(:,:,:,runnum) = all_data(m.y(1):m.y(2),m.x(1):m.x(2),3:3:end-1);
    %disp(['finished ' mat2str(runnum) ' out of ' mat2str(length(m.ttu)) ' runs.'])
    summer=summer+1;
    waitbar(runnum/numel(m.ttu),h,['finished ' mat2str(runnum) ' out of ' mat2str(length(m.ttu)) ' runs.']);
end
close(h)
data.b=rot90(blue_avg/summer,m.nrot);
data.g=rot90(green_avg/summer,m.nrot);
data.r=rot90(red_avg/summer,m.nrot);
datads.b=rot90(datads.b,m.nrot);
datads.g=rot90(datads.g,m.nrot);
datads.r=rot90(datads.r,m.nrot);
a = toc;
disp(['Total time - ' mat2str(a) ' seconds, taking about ' mat2str(a/numel(m.ttu)) ' seconds per trial'])

