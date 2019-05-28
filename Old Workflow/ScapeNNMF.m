%% Location is /local_mount/space/lfoivault/lfoivault1/SCAPE_DATA_BACKUP/clay_dendrites_gcamp6f_mouse1_040815/tiffs/runJ_stim_tiffs_diagshiftOW1
filepath = uigetdir('/local_mount/space/lfoivault/lfoivault1/SCAPE_DATA_BACKUP/clay_dendrites_gcamp6f_mouse1_040815/tiffs');
cd(filepath);

files = dir('gray*.tiff');
[~,ind]=sort({files.name});
files = files(ind);
%%
clear MIPdata
h = waitbar(0,'Loading images...');
for i = 1:395;
   imname = ['gray_runJ_stimfr' mat2str(i+5) '.tiff'];
   MIJ.run('Open...', ['path=['  strrep(filepath,'/','//') '//' imname ']']);
   %MIJ.run('Reslice [/]...', 'output=1.000 start=Left avoid');
   %MIJ.run('Open...', ['path=[//local_mount//space//lfoivault//lfoivault1//SCAPE_DATA_BACKUP//clay_dendrites_gcamp6f_mouse1_040815//tiffs//runJ_stim_tiffs_diagshiftOW1//' imname ']']);
   %MIJ.run('Z Project...', 'projection=[Average Intensity]');
   %MIPdata(:,:,i) = MIJ.getImage(['AVG_' imname]);
   data(:,:,:,i) = MIJ.getImage(imname);
   if i == 1
       %MIPdata = zeros(size(MIPdata,1),size(MIPdata,2),numel(files));
       data = zeros(size(data,1),size(data,2),size(data,3),numel(files));
       %MIPdata(:,:,i) = MIJ.getImage(['AVG_' imname]);
       data(:,:,:,i) = MIJ.getImage(imname);
   end
   
   MIJ.run('Close All')
   waitbar(i/numel(files),h);
end
close(h)

%% Making MIPS
% Top down
MIPtop = squeeze(mean(data,3));
MIPleft = squeeze(mean(data,2));
MIPside = squeeze(mean(data,1));

%% NNMF
MIPdata = MIPtop; %MIPdata = MIPdata(:,:,1:250);
ss(1) = size(MIPdata,1); ss(2) = size(MIPdata,2); ss(3) = size(MIPdata,3);
%W = Wt;
%[W,H] = nnmf(reshape(single(MIPdata),[ss(1)*ss(2) ss(3)]),24);
MIPrecon = reshape(W*H,size(MIPdata));
MIPresid = MIPrecon - MIPdata;
% For left MIP
%Hl = reshape(single(MIPdata),[ss(1)*ss(2) ss(3)])'/W';
%Wl = reshape(single(MIPdata),[ss(1)*ss(2) ss(3)])/H;
%%
figure
W = W/max(max(W));
cmap = hsv(size(H,1));
imgfull = zeros([ss(1) ss(2) 3]);
imgfulltop = zeros([ss(1) ss(2) 3]);
%imgfullleft = zeros([ss(1) ss(2) 3]);
for i = 1:24
    subplot(6,4,i)
    imgtmp = reshape(W(:,i)*cmap(i,:),[ss(1) ss(2) 3]);
    %imgtmp(imgtmp<.005) = 0;
    imagesc(rot90(imgtmp/max(W(:,i)),1)*10)
    title(mat2str(i))
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
end

%% Comp
% Make H non-negative
clear comp*
compstouse = [3 5 6 8 10 11 12 13 14 15 16 19 20 21 23 24];
ncomps = numel(compstouse);
comp = zeros([ss(1) ss(2) ss(3) ncomps],'single');

h = waitbar(0,'Making movie components...');
a = 1;
for i = compstouse
    comp(:,:,:,a) = reshape(W(:,i)*H(i,:),[ss(1) ss(2) ss(3)]);
    comp(:,:,:,a) = comp(:,:,:,a) 
    waitbar(i/ncomps,h)
    a = a+1;
end
close(h)

% Make NNMF movie
compvid = zeros([ss(1) ss(2) 3 ss(3)],'single');
h = waitbar(0,'Making movie frames...');
for i = 1:size(comp,3)
    compvid(:,:,:,i) = reshape(reshape(squeeze(comp(:,:,i,:)),[ss(1)*ss(2) ncomps])*hsv(ncomps),[ss(1) ss(2) 3]);
    waitbar(i/size(comp,3),h)
end
close(h)

clear M
brightness = 10/mean(max(max(max(compvid,[],1),[],2),[],4));
compvid_t = compvid;
compvid_t(compvid_t<.005) = 0;
%%
clear M
close all
% Top
%fig = figure('Color','Black','Position',[0 0 353*2 301*2]);
% Side
fig = figure('Color','Black','Position',[500 500 1200 400]);

for i = 1:ss(3)   
    imagesc(rot90(squeeze(compvid_t(:,:,:,i))*brightness,1)); daspect([2.9 2.2 1]);  
    axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    title(mat2str(i/10),'Color','w')
    M(i) = getframe(gcf);
end
close(fig)


%%
% Randy
%y_umPerPix = 1.384;
%z_umPerPix = 1.11;
%x_umPerPix = 313.51*info.daq.scanAngle/(info.daq.pixelsPerLine-2);

% Clay
y_umPerPix =  2.9;
z_umPerPix = 3.1;
x_umPerPix = 2.2;


%%

data = SCAPE_data(375-70:375+69,5:144,:,:);
clear SCAPE_data
for i = 1:size(data,2)
    data(:,i,:,:) = smooth3(squeeze(data(:,i,:,:)),'gaussian',[5 5 5]);
    disp(i)
end
ss = size(data);
bg = min(data(:,:,:,100:130),[],4);
data = data-repmat(bg,[1 1 1 ss(4)]);
[IDX,~] = kmeans(reshape(double(data),[ss(1)*ss(2)*ss(3),ss(4)]),18,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','correlation');







