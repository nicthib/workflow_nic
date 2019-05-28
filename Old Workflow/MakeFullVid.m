%% Initialize
addpath(genpath('/local_mount/space/enterprise/4/Personal Folders/Nic/'))
AFP = '/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/';
cd(AFP) % AFP is the Analysis Folder Path - where components, videos and MIDI files are stored.
% The m struct is used to store all run info, file paths, and analysis
% parameters. Basically, m stores all small variables that are neccesary to
% analyze runs.
m.camera = 'andor'; m.nrot = 3; % set m.firsttime to 1 if you haven't created metadata files.
                                                 % nrot is for rotating the images straight up.
m = GetMetaData(m.camera,0,0,0);
m.FileName = [m.mouse '_' m.run]; 
[m,data_raw] = LoadData(m); % First significant processing task

%% Masking
cd(AFP)
try % Tries to go into the mouse component dir - creates neccesary folders if it doesn't exist.
    cd(m.mouse)
catch
    mkdir(AFP,m.mouse); cd(m.mouse)
    mkdir([AFP m.mouse],'Videos')
    mkdir([AFP m.mouse],'MIDI')
    disp(['A folder was created for ' m.mouse]);
end

try % The params file contains the masks BW and BW_u - BW is the mask, and BW_u is the unilateral one for NNMF
    load([m.mouse '_params.mat'], 'BW','BW_u')
catch
    f = warndlg('No mask exists for this trial. You must create it before proceeding','Oops!');
    waitfor(f); clear f
    imagesc(mean(data.blue,3)); axis image
    h = impoly(gca); [XY1] = getPosition(h); X1 = XY1(:,1); Y1 = XY1(:,2);
    h = impoly(gca); [XY2] = getPosition(h); X2 = XY2(:,1); Y2 = XY2(:,2);
    BW = poly2mask(X1,Y1,m.height,m.width)+poly2mask(X2,Y2,m.height,m.width);
    BW(BW==0) = NaN;
    close all; clear h X* Y*; save([m.mouse '_mask.mat'],'BW');
end

%% Apply mask, convert, get seeds.
m.sz = size(data_raw.b,1); 
m.bs = 20/2^(9-log2(m.sz)); % This is entirely too complicated but I don't care.
m.interval = 1:(m.nFrames/m.nLEDs);
data_raw.b = data_raw.b.*repmat(BW,[1 1 size(data_raw.b,3)]);
data_raw.r = data_raw.r.*repmat(BW,[1 1 size(data_raw.b,3)]);
data_raw.g = data_raw.g.*repmat(BW,[1 1 size(data_raw.b,3)]);

[data.chbo,data.chbr,data.chbt] = convert_mariel(data_raw.g,data_raw.r,'g','r',m.interval,534);
data.gcamp=GcampMcCorrection(data_raw.b,data.chbr,data.chbo,m.interval,m.mua1,m.mua2);
clear data_raw data.chbo data.chbr
%%
cd([AFP m.mouse])
load([m.mouse '_params.mat'], 'seedpix')
if ~exist('seedpix')
    % Here, we create a downsampled version of the GCAMP data set and apply
    % NNMF to it in order to inform us of what seed regions to pick for
    % NNLS aalysis. 
    gcamp_tmp = reshape(imresize(data.gcamp(:,:,1:1000),[128 128]),[128^2 1000]);
    gcamp_tmp(isnan(gcamp_tmp)) = 0;
    [W_tmp,~] = nnmf(gcamp_tmp,12);
    [seedpix] = makeseedpix(W_tmp,12,0,0.2);
    seedpix = seedpix*m.sz/sqrt(size(gcamp_tmp,1)); % We can't use downsampled coordinates so we need to multiply them by the downsampling factor.
end
% LSQ analysis is here. We use a smooth subtraction to reduce background
% noise, and a smooth subtraction in the GCAMP signal to remove large
% drifts from the signal.
[H_gcamp_l,W_gcamp_l] = LSQanalysis(smooth3(data.gcamp,'box',[3 3 1])-smooth3(data.gcamp,'box',[3 3 101]),round(seedpix),m.bs);
[H_chbt_l,W_chbt_l] = LSQanalysis(smooth3(data.chbt,'box',[3 3 1])-smooth3(data.chbt,'box',[3 3 301]),round(seedpix),m.bs);
save([m.mouse '_' m.datatype],'H_*','W_*','m','seedpix');
%%
% Make Movies
m.compstouse = 1:12; ncomps = numel(m.compstouse);
cmap = hsv(ncomps); % cmap defines the colors you want to use for components. I typically use hsv.
cmap = cmap([1 18 2 17 3 16 4 15 5 14 6 13 7 12 8 11 9 10],:); %rearranging colors just as an experiment
for j = 1:2
    % H and W are the current ones being analyzed. Any changes I make to
    % them won't be kept on the original variables that are saved.
    if j == 1
        choice = 'gcamp_l';
        H = H_gcamp_l; W = W_gcamp_l;
    else
        choice = 'chbt_l';
        H = H_chbt_l; W = W_chbt_l;
    end
        
    H(H<0) = 0; % Remove negative values from H
    % comp represents each H and W component individually remixed. this 4D
    % matrix is neccesary for color remixing.
    comp = zeros([m.sz m.sz round(m.nFrames/m.nLEDs) ncomps],'single'); 
    h = waitbar(0,'Making movie components...'); a = 1;
    for i = m.compstouse
        comp(:,:,:,a) = reshape(W(:,i)*H(i,:),[m.sz m.sz round(m.nFrames/m.nLEDs)]);
        waitbar(i/ncomps,h); a = a+1;
    end; close(h)
    
    % compvid mixes the color coded components back into a 3D matrix.
    compvid = zeros([m.sz m.sz 3 round(m.nFrames/m.nLEDs)],'single');
    h = waitbar(0,'Making movie frames...');
    for i = 1:size(comp,3)
        compvid(:,:,:,i) = reshape(reshape(squeeze(comp(:,:,i,:)),[m.sz*m.sz ncomps])*cmap,[m.sz m.sz 3]);
        waitbar(i/size(comp,3),h)
    end; close(h)
    
    fig = figure('Color','Black','Position',[0,0,700,700]);
    % We define the brightness as the mean of the max values across each of the three color channels, divided by 3.
    brightness = 5/mean(max(max(max(compvid,[],1),[],2),[],4)); 
    % The video interpolates the frames down to a 10 fps video for
    % simplicity. i*m.fr/30 does this.
    for i = 1:1799
        imagesc(squeeze(compvid(:,:,:,ceil(i*m.fr/30)))*brightness)
        axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); colormap(hsv);
        M(i) = getframe(gcf);
    end; close(fig)
    
    cd(['/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/' m.mouse '/Videos'])
    Filename = [m.mouse '_' m.datatype '_' choice '.avi'];
    vid = VideoWriter(Filename,'Uncompressed AVI'); vid.FrameRate = 10;
    open(vid); writeVideo(vid,M); close(vid)
end

%% Webcam stimCCD side
close all
side = 'CCD'; % put CCD or stimCCD depending on what side you are doing.
% Webcam stimCCD side
if strcmp(side,'CCD')
    cd([(m.STIMpath{m.ttu}) '/' m.stimlist{m.ttu}(find(m.stimlist{m.ttu}~='_')) '_webcam']);
    fr_file = dir('*.txt'); fr_file = fr_file.name; FID = fopen(fr_file,'r');
    tfr1 = fscanf(FID, 'Frame Rate : %f'); fclose(FID);
    % tfr is the true frame rate, which is neccesary in order to sync
    % the webcam videos to the rest of the videos
    
    for i=0:length(dir)-4
        webcam_movie(:,:,i+1)=rgb2gray(imread([num2str(i) '.jpg']));
        if i == 0
            webcam_movie = zeros(size(webcam_movie,1),size(webcam_movie,2),length(dir)-3);
        end
    end
elseif strcmp(side,'stimCCD')
    stimnum = regexp(m.stimlist{m.ttu},'\d','Match');
    cd([(m.CCDpath{m.ttu}) '/' m.stimlist{m.ttu}(1:4) '_webcam' stimnum{1}]);
    fr_file = dir('*.txt'); fr_file = fr_file.name; FID = fopen(fr_file,'r');
    tfr2 = fscanf(FID, 'Frame Rate : %f');
    for i=0:length(dir)-4
        webcam_movie(:,:,i+1)=rgb2gray(imread([m.stimlist{m.ttu}(1:4) num2str(i) '.jpg']));
        if i == 0
            webcam_movie = zeros(size(webcam_movie,1),size(webcam_movie,2),length(dir)-3);
        end
    end
end

fig = figure('Color','Black','Position',[0,0,size(webcam_movie,1),size(webcam_movie,2)]); h = waitbar(0,['Making eye vid for ' m.datatype]);
for i = 1:5400
    imagesc(webcam_movie(:,:,round(i*tfr1/30)))
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); colormap gray;
    M(i) = getframe(gcf);
    waitbar(i/5400,h)
end
close(fig); close(h)

cd(['/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/' m.mouse '/Videos'])
Filename = [m.mouse '_' m.datatype '_webcam' side '.avi'];
vid = VideoWriter(Filename,'Uncompressed AVI');
vid.FrameRate = 30;
open(vid)
writeVideo(vid,M)
close(vid)
%%
MIDI_GCAMP(H_gcamp_l,m)
MIDI_CHBT(H_chbt_l,m)
