%% STARTUP
javaaddpath /usr/local/MATLAB/R2016a/java/mij.jar
javaaddpath /usr/local/MATLAB/R2016a/java/ij.jar
addpath(genpath('/local_mount/space/enterprise/4/Personal Folders/Nic/'))
addpath /local_mount/space/enterprise/4/NewFiji/Fiji.app/scripts

%%
clearvars -EXCEPT m BW*; close all; clc
m.camera = 'andor'; m.firsttime = 0; m.nrot = 3; 
m = GetMetaData(m); m.splitrun = 128;
m = FindStimRuns(m);
m.mouse = m.pathofmetadata(regexp(m.pathofmetadata,'[^/]+$'):end);
m.mua1 = .6; m.mua2 = .7;
m.ttu = 4;
m.bs = 20; m.sz = 512;
m.datatype = m.stimlist{m.ttu};
m.FileName = [m.mouse '_' m.datatype];
%% Main Analysis Code
m.datatype = m.stimlist{m.ttu};
BW = imresize(BW,[m.sz m.sz]); BW_u = imresize(BW_u,[m.sz m.sz]);
[m,data_raw] = LoadData(m);
%[m,data_raw] = LoadData_2LEDs(m);
m.interval = 1:(m.nFrames/m.nLEDs);
if ~exist('BW') % Make Mask
    imagesc(mean(data_raw.b,3));
    axis square
    h = impoly(gca);
    [XY1] = getPosition(h); X1 = XY1(:,1); Y1 = XY1(:,2);
    h = impoly(gca);
    [XY2] = getPosition(h); X2 = XY2(:,1); Y2 = XY2(:,2);
    BW = poly2mask(X1,Y1,m.height,m.width)+poly2mask(X2,Y2,m.height,m.width);
    BW_u = poly2mask(X1,Y1,m.height,m.width);
    close all;
    clear h XY1 XY2
end
if ndims(BW) == 2
    data_raw.b = data_raw.b.*repmat(BW,[1 1 size(data_raw.b,3)]);
    try data_raw.r = data_raw.r.*repmat(BW,[1 1 size(data_raw.b,3)]); catch; end
    data_raw.g = data_raw.g.*repmat(BW,[1 1 size(data_raw.b,3)]);
else
    data_raw.b = data_raw.b.*BW;
    try data_raw.r = data_raw.r.*BW; catch; end
    data_raw.g = data_raw.g.*BW;
end

tic
[data.chbo,data.chbr,data.chbt] = convert_mariel(data_raw.g,data_raw.r,'g','r',m.interval,534);
data.gcamp=GcampMcCorrection(data_raw.b,data.chbr,data.chbo,m.interval,m.mua1,m.mua2); gcamp = data.gcamp;
%[gcamp, green] = GCaMPcorr_green(data_raw,m);
%%
% Binning
gcamp = zeros(m.sz,m.sz,m.nFrames/m.nLEDs);
chbt = zeros(m.sz,m.sz,m.nFrames/m.nLEDs);
for i = 1:m.sz
    for j = 1:m.sz
        gcamp(i,j,:) = nanmean(nanmean(data.gcamp(i*2-1:i*2,j*2-1:j*2,:),2),1);
        chbt(i,j,:) = nanmean(nanmean(data.chbt(i*2-1:i*2,j*2-1:j*2,:),2),1);
    end
end
gcamp = imresize(data.gcamp,[256 256],'Method','bilinear');
%chbt = imresize(data.chbt,[256 256],'Method','bilinear');
gcamp = smooth3(data.gcamp-smooth3(data.gcamp,'box',[1 1 101]));
gcamp = smooth3(data.gcamp-smooth3(data.gcamp,'box',[3 3 101]),'box',[3 3 1]);

%chbt = smooth3(chbt,'box',[1 1 3]);

%% NNMF and LSQ

gcamp(isnan(gcamp)) = 0;
green(isnan(green)) = 0;
%chbt(isnan(chbt)) = 0;
%[W_gcamp_n,H_gcamp_n] = nnmf(reshape(gcamp,[m.sz^2 m.nFrames/m.nLEDs]),18);
%[W_chbt_n,H_chbt_n] = nnmf(reshape(data.chbt,[m.sz^2 m.nFrames/m.nLEDs]),18);

%[W_gcamp_n,H_gcamp_n] = nnmfanalysis(gcamp,18,BW_u);
%[W_chbt_n,H_chbt_n] = nnmfanalysis(chbtsm,18,round(imresize(BW_u,[256 256])));

%[seedpix] = makeseedpix(W_gcamp_n,18,0);

[H_gcamp_l,W_gcamp_l] = LSQanalysis(data.gcamp-smooth3(data.gcamp,'box',[1 1 301]),round(seedpix),m.bs);
[H_chbt_l,W_chbt_l] = LSQanalysis(data.chbt--smooth3(data.chbt,'box',[1 1 301]),round(seedpix),m.bs);

%cd(['/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/' m.mouse])
%save([m.mouse '_' m.stimlist{m.ttu} '_all'],'W*','H*','m','seedpix')

%% Residual analysis
choice = 'gcamp_l';
eval(['H = H_' choice '; W = W_' choice ';'])
[resid,HW_gcamp_corr_n] = getresid(gcamp,W,H,m.sz);

choice = 'chbt_l';
eval(['H = H_' choice '; W = W_' choice ';'])
[HW_gcamp_corr_l] = getresid(gcamp,W,H,m.sz);

%%
HW_gcamp_corr = HW_gcamp_corr_n;
close all
cjet = jet(200);
cjet(1,:) = 0;
figure('Color','Black','Position',[100 100 600 600]);
HW_gcamp_corr(HW_gcamp_corr==0) = NaN;
imagesc(HW_gcamp_corr)
title(sprintf('Correlation - LSQ vs GCaMP (average = %0.3f)',nanmean(nanmean(HW_gcamp_corr))),'Color','w')

colormap(cjet)
caxis([0 1])
c = colorbar; set(c,'Color','w','FontSize',10)
c.Label.String = 'Pearson''s Correlation';
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);

%% Plotting just components
close all
figure('Color','White','Position',[100 100 1000 1000]);
cmap = hsv(18);
imgfull = zeros([m.sz m.sz 3]);
imgfull2 = zeros([m.sz m.sz]);
%se = strel('cube',5);
%border = imresize(BW-imerode(BW,se),[m.sz m.sz]);
for i = 1:size(W,2)
    subplot(3,6,i)
    imgtmp = reshape(W(:,i)*cmap(i,:),[m.sz m.sz 3]);
    imgtmp = imgtmp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]));
   % imgtmp = imgtmp/nanmean(nanmean(nanmean(imgtmp)));
    imgfull = imgfull + imgtmp;
    imgfull2 = imgfull2 + reshape(W(:,i),[m.sz m.sz]);
    imagesc(imgtmp*5/max(max(max(imgtmp))))
    caxis([0 1])
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    hold on;
    title(mat2str(i),'Color','k','FontSize',9)
%     try
%         scatter(seedpix(i,1),seedpix(i,2),'w','Filled')
%     catch
%     end

end

% figure('Color','White','Position',[100 100 500 500]);
% imagesc(imgfull)
% axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
% hold on
% scatter(seedpix(:,1)/2,seedpix(:,2)/2,30,'w','Filled')
%hgexport(1,['Indiv_' choice])
%hgexport(2,['Full_' choice])

%% Comp
% Make H non-negative
H(H<0) = 0;
% for i = 1:18
%     H(i,:) = smooth(H(i,:))';
% end
clear comp*
compstouse = 1:12; %[2 5:15 18];
ncomps = numel(compstouse);
comp = zeros([m.sz m.sz floor(m.nFrames/m.nLEDs) ncomps],'single');

h = waitbar(0,'Making movie components...');
a = 1;
for i = compstouse
    comp(:,:,:,a) = reshape(W(:,i)*H(i,:),[m.sz m.sz floor(m.nFrames/m.nLEDs)]);
    waitbar(i/ncomps,h)
    a = a+1;
end
close(h)

% Make NNMF movie
compvid = zeros([m.sz m.sz 3 floor(m.nFrames/m.nLEDs)],'single');
h = waitbar(0,'Making movie frames...');
for i = 1:size(comp,3)
    compvid(:,:,:,i) = reshape(reshape(squeeze(comp(:,:,i,:)),[m.sz*m.sz ncomps])*hsv(ncomps),[m.sz m.sz 3]);
    waitbar(i/size(comp,3),h)
end
close(h)

%%
clear M
brightness = 5/mean(max(max(max(compvid,[],1),[],2),[],4));
%gcampsm = smooth3(gcamp,'box',[1 1 3]);
%residsm = smooth3(resid,'box',[1 1 3]);
%gcampsm(find(gcampsm==0)) = NaN;
%%
close all
%fig = figure('Color','Black','Position',[0,0,600,600]);
gcampsm = smooth3(data.gcamp,'box',[3 3 7])-smooth3(data.gcamp,'box',[3 3 101]);
chbtsm = smooth3(data.chbt,'box',[3 3 7])-smooth3(data.chbt,'box',[3 3 101]);

%%
fig = figure('Color','Black','Position',[0,0,200,200]);

jetc = hsv;
jetc(1,:) = 0;
for i = 1:1799
    %    imagesc(squeeze(imrotate(padarray(compvid(:,:,:,ceil(i*m.fr/30)),[0 50],'pre'),-3))*brightness)
    %                                           ^ for ketamine ^
    imagesc(squeeze(compvid(:,:,:,ceil(i*m.fr/30))*brightness))
    imagesc(data.gcamp(:,:,ceil(i*m.fr/30)))
    %imagesc(squeeze(compvid(:,:,ceil(i*m.fr/30))))
    
    %subplot(211)
   % imagesc(resid(:,:,ceil(i*m.fr/30)));     
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); colormap gray%colormap(jetc);
    %caxis([-1 1]*1e-1)
    caxis([-.05 .05])
    %subplot(212)
    %imagesc(chbtsm(:,:,i));     
    %axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); colormap jet;
    
    M(i) = getframe(gcf);
end
close(fig)
%%
cd(['/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/' m.mouse])
Filename = [m.mouse '_' m.datatype '_gcampraw.avi'];
myvideo = VideoWriter(Filename,'Uncompressed AVI');
myvideo.FrameRate = 10;
open(myvideo)
writeVideo(myvideo,M)
close(myvideo)

%%
Miji

%% PUPIL ANALYSIS
m = PupilAnalysis2(m,2); % 1 for stim side, 2 for camera side

%% Plot Pupil Trace
figure
%close all
t1 = 1/30:1/30:180;
m.Araw = m.A;
m.A1 = hampel(m.A,20);
m.A1(m.A1==0) = NaN;
m.A1(m.A1<300) = NaN;
m.A1(end,:) = NaN;
for i = 1:numel(m.ttu)
    patch(t1,m.A1(:,i),'Black','EdgeAlpha',1);
    hold on
end

%% Running Analysis
webcamloc = 2;
m.datatype = m.stimlist{m.ttu};
[TCout(:,a)] = getwebcamrunTC(m,webcamloc);

%%
clear data_raw
h = waitbar(0,'Downsampling...');
for i = 1:length(m.ttu)
    [datads.chbo(:,:,:,i),datads.chbr(:,:,:,i),datads.chbt(:,:,:,i)] = convert_mariel(data_rawds.g(:,:,:,i),data_rawds.r(:,:,:,i),'g','r',m.interval,534);
    datads.gcamp(:,:,:,i) = GcampMcCorrection(data_rawds.b(:,:,:,i),datads.chbr(:,:,:,i),datads.chbo(:,:,:,i),m.interval,m.mua1,m.mua2);
    waitbar(i/numel(m.ttu),h,['finished ' mat2str(i) ' out of ' mat2str(length(m.ttu)) ' runs.']);
end
close(h); a = toc;
disp(['Total time - ' round(mat2str(toc)) ' seconds, taking about ' round(mat2str(a/numel(m.ttu))) ' seconds per trial'])
%%
h = figure('Color','White','Position',[100 100 1000 1000]);

subplot(341)
imagesc(nanstd(data.chbt,[],3));
caxis([-12 12]*1e-6)
cmap = jet; cmap(1,:) = 1;
smap = jet(18);
c = [3 7 16 17];

colormap(cmap);
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); hold on;
for i = 1:4
    scatter(seedpix(i,1),seedpix(i,2),'Filled','MarkerFaceColor',smap(c(i),:),'MarkerEdgeColor','k')
end

%
subplot(342)
imagesc(nanstd(data.gcamp,[],3));
caxis([-.1 .1])
cmap = jet; cmap(1,:) = 1;
colormap(cmap);
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); hold on
for i = 1:4
    scatter(seedpix(i,1),seedpix(i,2),'Filled','MarkerFaceColor',smap(c(i),:),'MarkerEdgeColor','k')
end


imgfull = zeros([m.sz m.sz 3]);
imgfull2 = zeros([m.sz m.sz]);
for i = 1:4
    subplot(3,4,i+2)
    imgtmp = reshape(W_gcamp_l(:,i)*smap(c(i),:),[m.sz m.sz 3]);
    imgtmp = imgtmp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]));
    imagesc(imgtmp)
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);   
    
    
    subplot(3,4,i+2+4)
    imgtmp = reshape(W_chbt_l(:,i)*smap(c(i),:),[m.sz m.sz 3]);
    imgtmp = imgtmp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]));
    imagesc(imgtmp)
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);   

end


