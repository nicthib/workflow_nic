%% Only for first time
clear; close all; clc

[CCDpath,stimlist,pathofsave,STIMpath]=GetStimList(1);
camera='andor';
for stimcount=1:length(stimlist);
    warning off
    [saved]=SaveMetaData(camera,stimlist{stimcount},CCDpath{stimcount},pathofsave,1,STIMpath{stimcount});
    disp('Yes')
end

%% Set up paths again, select all runs in pathofmetadata
clear; close all; clc
addpath /local_mount/space/enterprise/4/'Personal Folders'/Mohammed/'Matlab Code'
addpath /local_mount/space/enterprise/4/'Personal Folders'/Mohammed/'Matlab Code'/AwakeAnalysis/
[CCDpath,stimlist,pathofmetadata,STIMpath]=GetStimList(1);
cd(pathofmetadata);
ccd_dir = dir;
ccd_dir = {ccd_dir.name};
subsetted=cellfun(@(x) x(strfind(x,'stim')),ccd_dir,'UniformOutput',false);
subsetted=cellfun(@isempty,subsetted,'UniformOutput',false);
subsetted=cell2mat(subsetted);
subsetted=ones(size(subsetted))-subsetted;
subsetted=logical(subsetted);
ccd_dir=(ccd_dir(subsetted));
clear subsetted
tiffdatflag=1;

metadataruns=ccd_dir;
%% index all stim runs
stim1_indices=[];
stim3_indices=[];
stim5_indices=[];
stim7_indices=[];
for i=1:length(CCDpath);
    i
    file_to_use=matfile([pathofmetadata '/' metadataruns{i}]);
    if file_to_use.tstim==5
        stim5_indices=[stim5_indices,i];
    elseif file_to_use.tstim==1
        stim1_indices=[stim1_indices,i];
    elseif file_to_use.tstim==3
        stim3_indices=[stim3_indices,i];
    elseif file_to_use.tstim==7
        stim7_indices=[stim7_indices,i];
    else
    end
end
%% index w/o stim data
stim5_indices=listdlg('PromptString','Stim 5','ListString',stimlist);
%%
stim1_indices=listdlg('PromptString','Stim 1','ListString',stimlist);

%% Check registration of images across all trials
height=512; width=512; tiffdatflag=1;
warning off; clear all_blue all_green all_red
for runnum=1:length(CCDpath)
    try
        runnum
        if tiffdatflag==0
            runtoload=[CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) '.tif'];
            infoimage=imfinfo(runtoload);
            TifLink=Tiff(runtoload,'r');
            TifLink.setDirectory(1);
            all_blue(:,:,runnum)=TifLink.read();
            TifLink.setDirectory(2);
            all_green(:,:,runnum)=TifLink.read();
            TifLink.close();
        else
            fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(0),2)]) num2str(0) '.dat'],'r','l');
            all_blue(:,:,runnum) = fread(fid,[height width],'uint16','l');
            fclose(fid); clear fid;
            fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(1),2)]) num2str(1) '.dat'],'r','l');
            all_green(:,:,runnum) = fread(fid,[height width],'uint16','l');
            fclose(fid); clear fid;
            fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(1),2)]) num2str(2) '.dat'],'r','l');
            all_red(:,:,runnum) = fread(fid,[height width],'uint16','l');
            fclose(fid); clear fid;
        end
    end
end
%% continued
close all;
%for i=[110:120]
for i=stim1_indices(1):stim1_indices(60)
    imagesc(imgradient(all_blue(:,:,(i)))-imgradient(all_blue(:,:,(stim1_indices(1))))); colormap gray; title(i);
    %imagesc(all_blue(:,:,i)); title(i)
    %waitforbuttonpress
    %pause(.01);
    drawnow
    
    %if i==106; waitforbuttonpress; c=caxis; else; caxis(c); pause(.1); end;
end
%% create gradient magnitudes for each of these trials in blue and green
i=1;

for j=1:length(CCDpath);
    %for j=1:length(ruList);
    [bluemag(:,:,i),~]=imgradient(all_blue(:,:,j));
    %[greenmag(:,:,i),~]=imgradient(all_green(:,:,j));
    i=i+1;
end
%% register images...?
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations=3000;
i=1
for j=[109 ]%1:length(CCDpath);
    j
    tform1 = imregtform(bluemag(100:350,100:350,j), bluemag(100:350,100:350,15), 'similarity', optimizer, metric);
    tformblue{i}=tform1;
    bluereg(:,:,i)=imwarp(all_blue(:,:,j), tformblue{i}, 'OutputView',imref2d([512 512]));
    
    %tform2 = imregtform(greenmag(100:412,100:412,i), greenmag(100:412,100:412,end), 'similarity', optimizer, metric);
    %tformgreen{i}=tform2;
    %greenreg(:,:,i)=imwarp(all_green(:,:,j), tformgreen{i}, 'OutputView',imref2d([512 512]));
    i=i+1;
    
end
%% load a resting state run
tic
runnum=60;
runtoload=[CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) '.tif'];
run('/local_mount/space/enterprise/4/Personal Folders/Mohammed/Matlab Code/AwakeAnalysis/LoadTiffRun.m')
toc
%% STIMULUS DATA: create green signal averages of your choice
tic
_use=stim5_indices(6:end); % YOUR CHOICE!
clear blue_pre_all green_pre_all red_pre_all
delete(gcp('nocreate'));
%parpool(8)
if tiffdatflag==0
    for i=1:length(trials_to_use)
        i
        runtoload=[CCDpath{trials_to_use(i)} '/' stimlist{trials_to_use(i)} '/' stimlist{trials_to_use(i)}(1:4) '.tif'];
        if i==1
            imginfo=imfinfo(runtoload);
            nFrames=length(imginfo);
            load(cell2mat(fullfile(pathofmetadata,metadataruns(trials_to_use(i)))))
            %frameRate=31.22; tpre=9.66; tstim=5; tpost=30-tpre-tstim;
            frame_pre=ceil(tpre*frameRate)+mod(floor(tpre*frameRate),3)+1; %DIRTYY
            frame_stim=floor((tpre+tstim)*frameRate);
        end
        blue_pre=imread(runtoload,'Index',frame_pre-6,'Info',imginfo);
        blue_pre_all(:,:,i)=blue_pre;
        green_pre=imread(runtoload,'Index',frame_pre-5,'Info',imginfo);
        green_pre_all(:,:,i)=green_pre;
        red_pre=imread(runtoload,'Index',frame_pre-4,'Info',imginfo);
        red_pre_all(:,:,i)=red_pre;
        parfor j=frame_pre:frame_stim
            all_data(:,:,j)=imread(runtoload,'Index',j,'Info',imginfo);
        end
        blue_diff(:,:,i)=double(mean(all_data(:,:,frame_pre:3:frame_pre+10),3))-double(blue_pre);
        green_diff(:,:,i)=double(mean(all_data(:,:,frame_pre+1:3:frame_stim),3))-double(green_pre);
        red_diff(:,:,i)=double(mean(all_data(:,:,frame_pre+2:3:frame_stim),3))-double(red_pre);
    end
    clear all_data
else
    for runnum=1:length(trials_to_use)  
        runnum
        runtoload=[CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/'];
        if runnum==1
            %imginfo=imfinfo(runtoload);
            %nFrames=length(imginfo);
            nFrames=length(dir(runtoload))-2;
            load(cell2mat(fullfile(pathofmetadata,metadataruns(trials_to_use(runnum)))))
            %frameRate=frameRate/3;
            %frameRate=31.22; tpre=9.66; tstim=5; tpost=30-tpre-tstim;
            frame_pre=ceil(tpre*frameRate)+mod(floor(tpre*frameRate),3); %DIRTYY
            frame_stim=floor((tpre+tstim)*frameRate);
            
            
        end
        %fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(0),2)]) num2str(0) '.dat'],'r','l');
        %imagesc(fread(fid,[height,width],'uint16','l')); waitforbuttonpress;
        fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(frame_pre-3),2)]) num2str(frame_pre-3) '.dat'],'r','l');
        blue_pre = fread(fid,[height width],'uint16','l');
        blue_pre_all(:,:,runnum)=blue_pre;
        fclose(fid); clear fid;
        fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(frame_pre-2),2)]) num2str(frame_pre-2) '.dat'],'r','l');
        green_pre = fread(fid,[height width],'uint16','l');
        green_pre_all(:,:,runnum)=green_pre;
        fclose(fid); clear fid;
        fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(frame_pre-1),2)]) num2str(frame_pre-1) '.dat'],'r','l');
        red_pre = fread(fid,[height width],'uint16','l');
        red_pre_all(:,:,runnum)=red_pre;
        fclose(fid); clear fid;
        parfor j=frame_pre-5:frame_stim+45
            fid = fopen([CCDpath{runnum} '/' stimlist{runnum} '/' stimlist{runnum}(1:4) repmat('0',[1,10-size(num2str(j),2)]) num2str(j) '.dat'],'r','l');
            all_data(:,:,j) = fread(fid,[height width],'uint16','l');
            fclose(fid); %clear fid;
        end
        blue_diff(:,:,runnum)=double(mean(all_data(:,:,frame_pre-3:3:frame_pre+20),3))-double(blue_pre);
        green_diff(:,:,runnum)=double(mean(all_data(:,:,frame_stim-1+20:3:frame_stim-1+45),3))-double(green_pre);
        red_diff(:,:,runnum)=double(mean(all_data(:,:,frame_pre-1:3:frame_stim),3))-double(red_pre);
    end
    clear all_data
end
delete(gcp('nocreate'))
toc
%% plot those avareges
blue_diff_cum=cumsum(blue_diff,3);%./cumsum(1:size(blue_diff,3));
green_diff_cum=cumsum(green_diff,3);%./cumsum(1:size(green_diff,3));
red_diff_cum=cumsum(red_diff,3);%./cumsum(1:size(red_diff,3));
%% visualize!!
close all;
figure;
for i=1:size(blue_diff,3)
    imagesc(blue_diff(:,:,i)); title(i); colormap jet;
    if i==1
        c=caxis; waitforbuttonpress;
    else
        caxis(c/5);
        pause(.1)
    end
end
%% STIMLUS: load and average all wanted trials
%trials_to_use=20;%([1:3 8:38 40 42:43 46:38 50:51 53:60]); % YOUR CHOICE!
%trials_to_use=stim5_indices(stim5_not_running_during_stim_percent(1:60)>.5);%stim5_indices(1:60);%stim5_indices_no_running(1:60);
%trials_to_use=stim5_indices(find(stim5_not_running_during_stim_percent>.4 & stim5_not_running_during_stim_percent2>.5));
%trials_to_use=4:8;
trials_to_use=stim5_indices(1:60);
%trials_to_use=;
%trials_to_use=trials_to_use_5s_pre;
%trials_to_use=(stim5_trials_to_keep);%
%trials_to_use=stim1_indices(57:end);
%trials_to_use=trials_to_keep_1s_post%stim5_indices(61:end);%trials_to_use_1s_post;%stim1_indices(61:end);%trials_to_use_5s;%trials_to_use_3s([1:9 11:end]);%stim5_indices(5:end);%stim5_indices(6:end);
%trials_to_use=9;
%trials_to_use=trials_to_use(stim1_trials_to_keep);
%trials_to_use=trials_to_use(even_less_trials);
%trials_to_use=cell2mat(stimlist(end-7));
%trials_to_use=trials_to_use([1:12 21:end]);
%%
runnum=1
if ~exist('NN')
    NN=1;
end
clear tmp_gcamp tmp_chbt tmp_chbo tmp_chbr
runtoload=[CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/']
nFrames=length(dir(runtoload))-2;
load(cell2mat(fullfile(pathofmetadata,metadataruns(trials_to_use(runnum)))))
nFrames=length(dir(runtoload))-2;
startrange=floor(tpre*(frameRate/3))-5:floor(tpre*(frameRate/3));
delete(gcp('nocreate'));
all_data=zeros(height,width,nFrames);
blue_avg=zeros(height,width,floor(nFrames/3));
green_avg=zeros(height,width,floor(nFrames/3));
red_avg=zeros(height,width,floor(nFrames/3));
%chbo_avg=zeros(height,width,floor(nFrames/3));
%chbr_avg=zeros(height,width,floor(nFrames/3));
%gcamp_avg=zeros(height,width,floor(nFrames/3));
%parpool(12)
%tic
summer=0;
if tiffdatflag==0
    %     for i=1:length(trials_to_use)
    %         i
    %         runtoload=[CCDpath{trial  s_to_use(i)} '/' stimlist{trials_to_use(i)} '/' stimlist{trials_to_use(i)}(1:4) '.tif'];
    %         if i==1
    %             frameRate=31.22; tpre=9.66; tstim=5; tpost=30-tpre-tstim;
    %             imginfo=imfinfo(runtoload);
    %             nFrames=length(imginfo);
    %             %load(cell2mat(fullfile(pathofmetadata,metadataruns(trials_to_use(i)))))
    %         end
    %         parfor j=1:nFrames
    %             all_data(:,:,j)=double(imread(runtoload,'Index',j,'Info',imginfo));
    %         end
    %         LED1=all_data(:,:,1:2:end); blue=blue(:,:,1:end-1);
    %         LED2=all_data(:,:,2::end);
    %         %LED3=all_data(:,:,3:3:end); [chbo,chbr,~] =
    %         %convert_mariel(green,red,'g','r',startrange,534);
    %         for k=1:1
    %             if 0
    %             addpath('/local_mount/space/juno/1/ConvertFiles/conversion_loadthesefiles');
    %             load dpffsMW_400to700nm.mat
    %             load Hb_spectra.mat
    %             addpath('/local_mount/space/juno/1/ConvertFiles/LED_spectra')
    %             load LED_spectra_0711
    %                 load('/local_mount/space/enterprise/4/Personal Folders/Teresa/m345_spectrum.mat')
    %                 %spectra_red=RED;
    %             % measured LED spectra:
    %             x530 = spectra_green(:,1);
    %             y530 = spectra_green(:,2);
    %             x470 = spectra_blue(:,1);
    %             y470 = spectra_blue(:,2);
    %             x630 = spectra_red(:,1);
    %             y630 = spectra_red(:,2);
    %             Lolim = 400;
    %             Hilim = 700;
    %
    %             splineyHb = spline([Lolim:2:Hilim],Hb(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
    %             splineyHbO = spline([Lolim:2:Hilim],Hb02(find(lambda==Lolim):find(lambda==Hilim)),[Lolim:0.5:Hilim]);
    %             spliney470 = spline(x470,y470,[Lolim:0.5:Hilim]);
    %             spliney530 = spline(x530,y530,[Lolim:0.5:Hilim]);
    %             spliney630 = spline(x630,y630,[Lolim:0.5:Hilim]);
    %
    %             splineydpff_488 = spline([Lolim:2:Hilim],dpff_488(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
    %             splineydpff_530 = spline([Lolim:2:Hilim],dpff_530(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
    %             splineydpff_630 = spline([Lolim:2:Hilim],dpff_630(find(waves==Lolim):find(waves==Hilim)),[Lolim:0.5:Hilim]);
    %
    %             % ~ area under curves for E ( do also for dpff) order is blue,
    %             % green, red
    %             EHb(1) = sum((1/sum(spliney470))*spliney470.*splineyHb);
    %             EHb(2) = sum((1/sum(spliney530))*spliney530.*splineyHb);
    %             EHb(3) = sum((1/sum(spliney630))*spliney630.*splineyHb);
    %             EHbO(1) = sum((1/sum(spliney470))*spliney470.*splineyHbO);
    %             EHbO(2) = sum((1/sum(spliney530))*spliney530.*splineyHbO);
    %             EHbO(3) = sum((1/sum(spliney630))*spliney630.*splineyHbO);
    %
    %             %still need to incorporate this into the DPF
    %             DPF(1) = sum((1/sum(spliney470))*spliney470.*splineydpff_488);
    %             DPF(2) = sum((1/sum(spliney530))*spliney530.*splineydpff_530);
    %             DPF(3) = sum((1/sum(spliney630))*spliney630.*splineydpff_630);
    %
    %             mua_b(:,:,:)= EHb(1)*chbr+EHbO(1)*chbo;
    %             mua_g(:,:,:) = EHb(2)*chbr+EHbO(2)*chbo;
    %             blue_adjusted = exp(-DPF(1)*mua_b(:,:,:));
    %             clear mua_b mua_g
    %             gcamp = (blue./repmat(mean(blue(:,:, startrange),3),[1,1,size(blue,3)]))./sqrt(blue_adjusted)./sqrt(green./repmat(mean(green(:,:, startrange),3),[1,1,size(green,3)]));
    %             gcamp = gcamp-1;
    %             end
    %         end
    %         %blue_avg=blue_avg+blue; green_avg=green_avg+green;
    %         %red_avg=red_avg+red; chbo_avg=chbo_avg+chbo;
    %         %chbr_avg=chbr_avg+chbr; gcamp_avg=gcamp_avg+gcamp;
    %     end
else
    for runnum=1:length(trials_to_use)
        %runnum
        %frameRate=31.22; tpre=9.66; tstim=5; tpost=30-tpre-tstim;
        %if runnum==1
        %runtoload=[CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/'];
        %nFrames=length(dir(runtoload))-2;
        %load(cell2mat(fullfile(pathofmetadata,metadataruns(trials_to_use(runnum)))))
        %frameRate=31.22; tpre=9.66; tstim=5; tpost=30-tpre-tstim;
        %end
        disp([CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)}(1:4)]);
        %Par here
        for j=0:nFrames-1
            %    j
            %disp([CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)}(1:4) repmat('0',[1,10-size(num2str(j),2)]) num2str(j) '.dat']);
            %if runnum<37
                fid = fopen([CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)}(1:4) repmat('0',[1,10-size(num2str(j),2)]) num2str(j) '.dat'],'r','l');
            %else
            %    fid = fopen([CCDpath{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)} '/' stimlist{trials_to_use(runnum)}(1:4) '1' repmat('0',[1,10-size(num2str(j),2)]) num2str(j) '.dat'],'r','l');
            %end
            all_data(:,:,j+1) = fread(fid,[height width],'uint16','l');
            fclose(fid); %clear fid;
        end
        %blue=all_data(:,:,1:3:end); blue=blue(:,:,1:end-1);
        %green=all_data(:,:,2:3:end);
        %red=all_data(:,:,3:3:end);
        %[chbo,chbr,~] = convert_mariel(green,red,'g','r',startrange,534);
        %NN=1;
        
        if 0 %trials_to_use(runnum)>=6 && trials_to_use(runnum)<18
            
            blue_avg=blue_avg+imwarp(all_data(:,:,1:3:end-1),tformblue{1},'OutputView',imref2d([512 512]));
            green_avg=green_avg+imwarp(all_data(:,:,2:3:end),tformblue{1},'OutputView',imref2d([512 512]));
            red_avg=red_avg+imwarp(all_data(:,:,3:3:end),tformblue{1},'OutputView',imref2d([512 512]));
        elseif 0 %trials_to_use(runnum)>=18 && trials_to_use(runnum)<133
            blue_avg=blue_avg+imwarp(all_data(:,:,1:3:end-1),tformblue{2},'OutputView',imref2d([512 512]));
            green_avg=green_avg+imwarp(all_data(:,:,2:3:end),tformblue{2},'OutputView',imref2d([512 512]));
            red_avg=red_avg+imwarp(all_data(:,:,3:3:end),tformblue{2},'OutputView',imref2d([512 512]));
        elseif 0 %trials_to_use(runnum)>=19 && trials_to_use(runnum)<133
            blue_avg=blue_avg+imwarp(all_data(:,:,1:3:end-1),tformblue{3},'OutputView',imref2d([512 512]));
            green_avg=green_avg+imwarp(all_data(:,:,2:3:end),tformblue{3},'OutputView',imref2d([512 512]));
            red_avg=red_avg+imwarp(all_data(:,:,3:3:end),tformblue{3},'OutputView',imref2d([512 512]));
        else
            10
            blue_avg=blue_avg+all_data(:,:,1:3:end-1);
            green_avg=green_avg+all_data(:,:,2:3:end);
            red_avg=red_avg+all_data(:,:,3:3:end);
            %blue_avg=all_data(:,:,1:3:end-1);
            %green_avg=all_data(:,:,2:3:end-1);
            %red_avg=all_data(:,:,3:3:end);
        end
        % 1
        for k=1:NN
            %trials_to_use(runnum)
            if 0 % trials_to_use(runnum)<112;
                %[chbo,chbr,chbt] = convert_mariel(all_data(y(1,k):y(1,k)+roi_size,x(1,k):x(1,k)+roi_si  ze,2:3:end),all_data(y(1,k):y(1,k)+roi_size,x(1,k):x(1,k)+roi_size,3:3:end),'g','r',startrange,534);
                %gcamp=GcampMcCorrection(all_data(y(1,k):y(1,k)+roi_size,x(1,k):x(1,k)+roi_size,1:3:end-1),chbr,chbo,startrange,.5,.5);
                %tmp_gcamp(:,runnum,k)=squeeze(mean(mean(gcamp(:,:,1:end))));
                %tmp_chbt(:,runnum,k)=squeeze(mean(mean(chbt(:,:,1:end))));
                %tmp_chbo(:,runnum,k)=squeeze(mean(mean(chbo(:,:,1:end))));
                %   tmp_chbr(:,runnum,k)=squeeze(mean(mean(chbr(:,:,1:end))));
                tmp_blueavg(:,trials_to_use(runnum),k)=squeeze(mean(mean(all_data(y1(1,k):y1(1,k)+roi_size,x1(1,k):x1(1,k)+roi_size,1:3:end-1))));
                 tmp_greenavg(:,trials_to_use(runnum),k)=squeeze(mean(mean(all_data(y1(1,k):y1(1,k)+roi_size,x1(1,k):x1(1,k)+roi_size,2:3:end))));
                 tmp_redavg(:,trials_to_use(runnum),k)=squeeze(mean(mean(all_data(y1(1,k):y1(1,k)+roi_size,x1(1,k):x1(1,k)+roi_size,3:3:end))));
            elseif 0
                 tmp_blueavg(:,(runnum),k)=squeeze(mean(mean(all_data(y(1,k):y(1,k)+roi_size,x(1,k):x(1,k)+roi_size,1:3:end-1))));
                 tmp_greenavg(:,(runnum),k)=squeeze(mean(mean(all_data(y(1,k):y(1,k)+roi_size,x(1,k):x(1,k)+roi_size,2:3:end))));
                tmp_redavg(:,(runnum),k)=squeeze(mean(mean(all_data(y(1,k):y(1,k)+roi_size,x(1,k):x(1,k)+roi_size,3:3:end))));
            end
        end
        %chbo_avg=chbo_avg+chbo;
        %chbr_avg=chbr_avg+chbr;
        %gcamp_avg=gcamp_avg+gcamp;
        %end
        %clear all_data;    
        %blue=all_data;
        %save(['cm42_6_' stimst{trials_to_use(runnum)} '_processed.mat'],'blue','green','-v7.3');
        %[chbo,chbr,chbt] = convert_mariel(green_avg,red_avg,'g','r',startrange,534);
        %im_b=blue_avg(:,:,100);
        %im_g=green_avg(:,:,100);
        %im_r=red_avg(:,:,100);
        %frameRate=(31.2207/3);
        %save(['new_cm42_4_' stimlist{trials_to_use(runnum)} '_processed.mat'],'chbt','frameRate','im_g','im_b','im_r','-v7.3');
        summer=summer+1;
    end
    
end
%clear all_data;
%toc
%delete(gcp('nocreate'));
i=summer;
%chbt_avg=chbo_avg+chbr_avg;
%chbo_avg=chbo_avg/i;
%chbr_avg=chbr_avg/i;
%chbt_avg=chbt_avg/i;
blue_avg=(blue_avg/i);
green_avg=(green_avg/i);
red_avg=(red_avg/i);
%blue_avg=(blue_avg);
%green_avg=(green_avg);
%red_avg=(red_avg);


%}
%gcamp_avg=gcamp_avg/i;
%clear i runnum

%% STIMULUS: visualize
close all;
for i=90:size(chbt_avg_pre,3)
    imagesc(cat(2,(gcamp_avg_pre(:,:,i))*.0005,chbt_avg_pre(:,:,i))); title(i*3/frameRate); colormap jet; caxis([-1 1]*1e-5); %colorbar
    %imagesc(gcamp_avg(:,:,i))%-blue_avg(:,:,89));
    %if i==90; c=caxis; else; caxis(c); end;
    %caxis([-1 1]*1e-1);
    pause(1*3/frameRate);
end
%%
[chbo_avg1,chbr_avg1,chbt_avg1] = convert_mariel(green_avg,red_avg,'g','r',95:100,534);
gcamp1=GcampMcCorrection(blue_avg,chbr_avg1,chbo_avg1,95:100,.5,.6);