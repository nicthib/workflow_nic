load('/local_mount/space/revault/revault1/Sam_SCAPE1p/mouse17_run3_1_forliam.mat')
bg = min(SCAPE_data(:,:,:,100:130),[],4);
ss = size(SCAPE_data);
SCAPE_data_bgs = SCAPE_data-repmat(min(SCAPE_data(:,:,:,100:130),[],4),[1 1 1 ss(4)]); clear SCAPE_data
SCAPE_data_bgs = reshape(SCAPE_data_bgs,[prod(ss(1:3)) ss(4)]);
SCAPE_data_bgs_ds = squeeze(mean(mean(mean(reshape(SCAPE_data_bgs(1:2*floor(ss(1)/2),1:2*floor(ss(2)/2),1:2*floor(ss(3)/2),:),[2,floor(ss(1)/2),2,floor(ss(2)/2),2,floor(ss(3)/2),ss(4)]),5),3),1));
ss2 = size(SCAPE_data_bgs_ds);
%%
IPmask = findimportantpix(reshape(SCAPE_data_bgs_ds,[prod(ss2(1:3)),ss2(4)]));
IPmask = reshape(IPmask,ss2(1:3)); IPmask(IPmask == 0) = NaN; 
[IDX,~] = kmeans(reshape(SCAPE_data_bgs_ds.*repmat(IPmask,[1 1 1 ss2(4)]),[prod(ss2(1:3)) ss(4)]),85,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','correlation');
IDX = reshape(IDX,ss2(1:3));
%% Visualize
for i = 1:max(IDX(:))
    [x,y,z] = ind2sub(ss2(1:3),find(IDX == i));
    shp = alphaShape(x,y,z,1);
    plot(shp,'EdgeColor','none','FaceColor',rand(1,3))
    hold on
end
%%
for i = 1:max(IDX(:))
    [x,y,z] = ind2sub(ss2(1:3),find(IDX == i));
    tmpH = zeros(ss2(4),1);
    for j = 1:numel(x)
        tmpH = tmpH + squeeze(squeeze(squeeze(mean(mean(mean(SCAPE_data_bgs(x(j)*2-1:x(j)*2,y(j)*2-1:y(j)*2,z(j)*2-1:z(j)*2,:),3),2),1))));
    end
    H(i,:) = tmpH/numel(x);
end
clear x y z tmpH shp
%% NNLS
SCAPE_data_bgs = reshape(SCAPE_data_bgs,[prod(ss(1:3)) ss(4)]);
W = zeros(prod(ss(1:3)),size(H,1));
tic
parfor i = 1:prod(ss(1:3)) % Perform NNLS. parfor makes this less hellish
    W(i,:) = lsqnonneg(H',squeeze(SCAPE_data_bgs(i,:))');    
end
toc
%% Sort components spatially
for i = 1:size(W,2)
    I = reshape(W(:,i),ss(1:3));
    frame = regionprops(true(size(I)),I, 'WeightedCentroid');
    cm(i,:) = frame.WeightedCentroid;
    i
end
[~,s] = sort(cm(:,2));
H_s = H(s,:);
W_s = W(:,s);

%% Visualize and choose components
hFig = figure('Position',[200 200 1000 700]);
cmap = hsv(size(H_s,1));
compstouse = [];
for i = 1:85
    subplot(221)
    SCAPE_small = SCAPE_data(1:ss(1)/2,:,1:ss(3)/2,:);
    imagesc(rot90(squeeze(max(reshape((W_s(:,i))*cmap(i,:),[ss(1:3) 3]),[],3)),1))
    axis image
    subplot(222)
    imagesc(rot90(squeeze(max(reshape((W_s(:,i))*cmap(i,:),[ss(1:3) 3]),[],1)),2))
    axis image
    subplot(223)
    imagesc(rot90(squeeze(max(reshape((W_s(:,i))*cmap(i,:),[ss(1:3) 3]),[],2)),3))
    axis image
    subplot(224)
    plot(H_s(i,:)./pbleach)
    title(['Component ' mat2str(i)])
    colormap gray

    a = waitforbuttonpress;
    if (a == 0)
        compstouse = [compstouse i];
    end
end
close all
H_final = H_s(compstouse,:)./repmat(pbleach,[numel(compstouse),1])-1;
W_final = W_s(:,compstouse);
%% MUSIC

load('/local_mount/space/revault/revault1/Sam_SCAPE1p/mouse17_run3_1_forliam.mat', 'info')
m.framerate = info.daq.scanRate;
m.keys = [1 3 5 6 8 10 12]; 
m.keys = [m.keys m.keys+12];
m.keys = [m.keys m.keys+24];
m.keys = [m.keys m.keys+48];
m.keys = m.keys+20;
m.nLEDs = 1;
MIDI_GCAMP(H_final,H_final,m,.3)
clear stack
%% Create stacks and save into tiffs.
stack = W_final*H_final;
smin = min(stack(:));
smax = max(stack(:));
clear stack
cmap = hsv(size(H_final,1));
h = waitbar(0,'Making stack...');
ss = [750 149 140 564];
for i = 1:ss(4)
    % load individal timepoints and convert to rgb tiff matrix.
    stack = zeros([ss(1:3) 3]);
    for j = 1:size(H_final,1)
        stack = stack+reshape(reshape(W_final(:,j)*H_final(j,i),[prod(ss(1:3)) 1])*cmap(j,:),[ss(1:3) 3]);
    end
    stack = stack-smin; stack = stack/smax;
    stack = uint16(round(stack*2^16));
    imgFileName = ['mouse17_run3_1_' mat2str(i) '.tiff'];
    for j = 1:size(stack,3)
        imwrite(squeeze(stack(:,:,j,:)), imgFileName, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
    end
    waitbar(i/ss(4),h)
end
close(h)

%% Make Stitched video
xcal = info.GUIcalFactors.x_umPerPix;
ycal = info.GUIcalFactors.y_umPerPix;
zcal = info.GUIcalFactors.z_umPerPix;
ss_i = [1033 175 351 564];
ss = [750 149 140 564];
cmap = hsv(numel(compstouse));
vidObj = VideoWriter('outputfull.avi','Uncompressed AVI');
vidObj.FrameRate = m.framerate;
open(vidObj)
fig = figure('Position',[0 0 ss_i(1)+ss_i(2) (ss_i(3)+ss_i(2))],'Color','k');
for i = 1:ss(4)
    frame = reshape(W_final*bsxfun(@times,cmap,H_final(:,i)/max(H_final(:))),[ss(1:3) 3]);
    frame = frame*2;
    for k = 1:3
        frame2(:,:,:,k) = interp3([0:ss(2)-1]*zcal,[0:ss(1)-1]*ycal,[0:ss(3)-1]'*xcal,squeeze(frame(:,:,:,k)),[0:ss(2)*zcal],[0:ss(1)*ycal],[0:ss(3)*xcal]','nearest');
    end
    figure(fig)
    im = [rot90(squeeze(max(frame2,[],2))) zeros(ss_i(3),ss_i(3),3) ; rot90(squeeze(max(frame2,[],3))) rot90(squeeze(max(frame2,[],1)),2)];
%     for k = 1:3
%         im(:,:,k) = medfilt2(im(:,:,k),[3 3]);
%     end
    imagesc(im)
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    writeVideo(vidObj,getframe(fig));
end
close(vidObj)




