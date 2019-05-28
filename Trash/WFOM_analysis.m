clear
addpath(genpath('/local_mount/space/juno/1/Software/MIAO/MIAO_v2'))
addpath('/local_mount/space/enterprise/4/Personal Folders/Nic/New Workflow')
%m.mouse = 'cm62_2'; 
m.mouse = 'cm72_9'; 
m.run = 'runG'; 
m.stim = 1;
m.dsf = 1;
%%
m.savepath = ['/local_mount/space/revault/revault2/MusicPaper_nic/Behavioral'];
m.ccddir = findmousefolder(m.mouse);
[m,data] = LoadDataPCA(fullfile(m.ccddir,m.run,[m.run '_stim_' mat2str(m.stim)]),'bodn',m.dsf,[.5 .6],m,100:200);
data.chbt = data.chbo + data.chbr; %data = rmfield(data,{'chbo','chbr'});
m.sz = size(data.gcamp,1);
data.gcamp = rot90(data.gcamp,2);
data.chbt = rot90(data.chbt,2);
% Preprocessing
load(fullfile(m.savepath,[m.mouse '_mask.mat']))
fprintf('Applying Mask...\n')
data.gcamp = data.gcamp.*repmat(round(imresize(BW,1/m.dsf)),[1 1 size(data.gcamp,3)]);
data.chbt = data.chbt.*repmat(round(imresize(BW,1/m.dsf)),[1 1 size(data.chbt,3)]);
disp('Smoothing gcamp...')
data.gcampsm = smooth3(data.gcamp,'box',[3 3 3])-smooth3(data.gcamp,'box',[3 3 101]);
disp('Smoothing chbt...')
data.chbtsm = smooth3(data.chbt,'box',[3 3 3]);
%%
[IDX,~] = kmeans(reshape(data.gcampsm,[size(data.gcampsm,1)*size(data.gcampsm,2) size(data.gcampsm,3)]),18,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','correlation');

H_gcamp = getHfromKmeans(data.gcampsm,IDX,BW,BW_u);
H_chbt = getHfromKmeans(data.chbtsm,IDX,BW,BW_u);

[~,W_gcamp] = LSQanalysis(data.gcampsm,H_gcamp);
[~,W_chbt] = LSQanalysis(data.chbtsm,H_chbt);

[H_gcamp,W_gcamp,H_chbt,W_chbt] = arrangecomps(H_gcamp,W_gcamp,H_chbt,W_chbt);


%%
for i = 1:18
    subplot(211)
    imagesc(reshape(W_gcamp(:,i),[512 512]))
    axis image
    subplot(212)
    imagesc(reshape(W_chbt(:,i),[512 512]))
    axis image
    waitforbuttonpress()
end
%% Make videos
m.compstokill = 7;
cmap = hsv(size(H_gcamp,1));
%cmap = cmap(randperm(numel(m.compstouse)),:);
m.interp = 1; % set to 1 if you want to interpolate to camera framerate.
m.brightness = 50;
MakeCompVid_new(H_gcamp,W_gcamp,m,'gcamp',cmap)
m.compstokill = [7 18];
m.brightness = 200000;
MakeCompVid_new(H_chbt,W_chbt, m,'chbt',cmap)

%% Make MIDI files
m.keys = [0 3 5 7 10]+36; m.keys = [m.keys m.keys+12 m.keys+24 m.keys+36 m.keys+48];
m.compstokill = 7;
MIDI_GCAMP(H_gcamp(m.compstouse,:),H_gcamp(m.compstouse,:),m,0.001)
m.compstokill = [7 18];
MIDI_CHBT(H_chbt,m)







