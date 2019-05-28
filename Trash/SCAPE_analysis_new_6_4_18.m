
load('/local_mount/space/revault/revault1/Sam_SCAPE1p/mouse17_run3_1_forliam.mat')
bg = min(SCAPE_data(:,:,:,100:130),[],4);
ss = size(SCAPE_data);
SCAPE_data_bgs = SCAPE_data-repmat(bg,[1 1 1 ss(4)]);
SCAPE_sm = zeros(ss);
parfor i = 1:ss(2)
    SCAPE_sm(:,i,:,:) = smooth3(squeeze(SCAPE_data_bgs(:,i,:,:)),'gaussian',[5 5 5]);
end
%% Smooth in time
bg = min(SCAPE_data(:,:,:,100:130),[],4);
ss = size(SCAPE_data);
SCAPE_data_bgs = SCAPE_data-repmat(bg,[1 1 1 ss(4)]);
clear SCAPE_data
SCAPE_data_bgs = reshape(SCAPE_data_bgs,[prod(ss(1:3)) ss(4)]);
SCAPE_data_bgs_s = zeros([prod(ss(1:3)) ss(4)]);
parfor i = 1:size(SCAPE_data_bgs,1)
    SCAPE_data_bgs_s(i,:) = smooth(double(SCAPE_data_bgs(i,:)),3);
end
SCAPE_data_bgs_s = reshape(SCAPE_data_bgs_s,ss);
SCAPE_data_bgs_s_ds = squeeze(mean(mean(mean(reshape(SCAPE_data_bgs_s(1:2*floor(ss(1)/2),1:2*floor(ss(2)/2),1:2*floor(ss(3)/2),:),[2,floor(ss(1)/2),2,floor(ss(2)/2),2,floor(ss(3)/2),ss(4)]),5),3),1));

%% Optional
SCAPE_sm = SCAPEsm(1:end,1:148,1:end,:);
% Downsample data for getting seed regions
SCAPE_sm2 = squeeze(mean(mean(mean(reshape(SCAPE_sm(1:2*floor(ss(1)/2),1:2*floor(ss(2)/2),1:2*floor(ss(3)/2),:),[2,floor(ss(1)/2),2,floor(ss(2)/2),2,floor(ss(3)/2),ss(4)]),5),3),1));
clear SCAPE_data SCAPE_data_bgs
%% NEW METHOD
% Take downsampled data, and for each voxel, smooth subtract and determine
% how many values rise above 3*SD of that TC. If enough do, it is
% considered an "interestign" voxel and is kept in a bianry mask. Then,
% k-means the binary mask into euclidean-seperated chunks, obtaining H from
% the full resolution data using these chunks.
P = findimportantpix(reshape(SCAPE_data_bgs_s_ds,[prod(ss2(1:3)),ss2(4)]));
P = reshape(P,ss2(1:3));
[a,b,c] = ind2sub(ss2(1:3),find(P));
[~,cnts] = kmeans([a b c],101,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','sqeuclidean');
[~,s] = sort(cnts(:,1));
[IDX,cnts] = kmeans(reshape(SCAPE_data_bgs_s_ds.*repmat(P,[1 1 1 ss2(4)]),[prod(ss2(1:3)) ss(4)]),85,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','correlation');

% IDX has seed #
P_new = P;
P_new(find(P)) = IDX;
%%
figure
for i = 1:max(IPmask(:))
    frame = IPmask;
    frame(frame ~= i) = 0;
    [a1,b1,c1] = ind2sub(ss2(1:3),find(frame));
    shp = alphaShape(a1,b1,c1,1);
    plot(shp,'EdgeColor','none','FaceColor',rand(1,3))
    hold on
    J{i} = frame/i;
end
%%
% J contains each seed region for getting TC
% Upscale J
clear H
bkgH = squeeze(squeeze(squeeze(mean(mean(mean(SCAPE_data_bgs_s,3),2),1))));
for i = 1:numel(J)
    frame = zeros([ss2(1:3)]);
    [x,y,z] = ind2sub(size(frame),find(J{i}));
    tmpH = zeros(ss2(4),1);
    for j = 1:numel(x)
        tmpH = tmpH+squeeze(squeeze(squeeze(mean(mean(mean(SCAPE_data_bgs_s(x(j)*2-1:x(j)*2,y(j)*2-1:y(j)*2,z(j)*2-1:z(j)*2,:),3),2),1))));
        %tmpH = tmpH+squeeze(squeeze(squeeze(mean(mean(mean(SCAPE_sm2(x(j),y(j),z(j),:),3),2),1))));
    end
    H(i,:) = tmpH/numel(x);
    i
end
% ADD BACKGROUND COMPONENT
H = [bkgH'; H];
%%
mail = 'dnt2111@gmail.com'; %Your GMail email address
password = 'G65414l@';  %Your GMail password
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
% Send the email.  Note that the first input is the address you are sending the email to
%% NNLS
SCAPE_data_bgs_s = reshape(SCAPE_data_bgs_s,[prod(ss(1:3)) ss(4)]);
W = zeros(prod(ss(1:3)),size(H,1));
parfor i = 1:prod(ss(1:3)) % Perform NNLS. parfor makes this less hellish
    W(i,:) = lsqnonneg(H',squeeze(SCAPE_data_bgs_s(i,:))');    
end
sendmail(mail,'LSQnonneg is DONE!','yaaaaaaaaaaaaaaaaaaaaaaaaaaaaay.')

%% Sort
for i = 1:size(W,2)
    I = reshape(W(:,i),ss(1:3));
    frame = regionprops(true(size(I)),I, 'WeightedCentroid');
    cm(i,:) = frame.WeightedCentroid;
    i
end
[~,s] = sort(cm(:,2));
H_s = H(s,:);
W_s = W(:,s);
%% Visualize
hFig = figure('Position',[200 200 1000 700]);
cmap = hsv(size(H_s,1));
compstouse = [];
for i = 1:85
    subplot(221)
    imagesc(rot90(squeeze(max(reshape((W_s(:,i))*cmap(i,:),[ss(1:3) 3]),[],3)),1))
    axis image
    subplot(222)
    imagesc(rot90(squeeze(max(reshape((W_s(:,i))*cmap(i,:),[ss(1:3) 3]),[],1)),2))
    axis image
    subplot(223)
    imagesc(rot90(squeeze(max(reshape((W_s(:,i))*cmap(i,:),[ss(1:3) 3]),[],2)),3))
    axis image
    subplot(224)
    plot(H_s(i,:))
    title(['Component ' mat2str(i)])
    colormap gray

    a = waitforbuttonpress;
    if (a == 0)
        compstouse = [compstouse i];
    end
end
%% MUSIC
load('/local_mount/space/revault/revault1/Sam_SCAPE1p/mouse17_run3_1_forliam.mat', 'info')
% Subtract background

m.framerate = info.daq.scanRate;
m.keys = [1 3 5 6 8 10 12]; 
m.keys = [m.keys m.keys+12];
m.keys = [m.keys m.keys+24];
m.keys = [m.keys m.keys+48];
m.keys = m.keys+20;

m.nLEDs = 1;

MIDI_GCAMP(H_final,H_final,m,.4)

%% Top and side movies
% Z - 149, Y - 750, X - 140
close all
vidObj = VideoWriter('SCAPE_daspect2.avi','Uncompressed AVI'); 
vidObj.FrameRate = m.framerate;
open(vidObj)

fig1 = figure('Position',[100 100 800 578]);
fig2 = figure('Position',[100 100 1600 800]);
fig3 = figure('Position',[100 100 1600 800]);
figfull = figure('Position',[0 0 1855 832]);

xcal = info.GUIcalFactors.x_umPerPix;
ycal = info.GUIcalFactors.y_umPerPix;
zcal = info.GUIcalFactors.z_umPerPix;
%W_final = W_s(:,compstouse);
%H_final = H_s(compstouse,:)-repmat(H_s(1,:),[numel(compstouse) 1]);
%H_final = H_final/max(H_final(:));

clear frame2
ss = [750 149 140 564];
cmap = hsv(size(H_final,1));

for i = 1:ss(4)
    frame = reshape(W_final*bsxfun(@times,cmap,H_final(:,i)),[ss(1:3) 3]);
    for k = 1:3
        frame2(:,:,:,k) = interp3([0:ss(2)-1]*zcal,[0:ss(1)-1]*ycal,[0:ss(3)-1]'*xcal,squeeze(frame(:,:,:,k)),[0:ss(2)*zcal],[0:ss(1)*ycal],[0:ss(3)*xcal]');
    end
    figure(fig1)
    imagesc(rot90(squeeze(max(frame2,[],1)),2))
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    M1 = getframe(gca);
    
    figure(fig2)
    imagesc(rot90(squeeze(max(frame2,[],2)),1))
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    M2 = getframe(gca);
    
    figure(fig3)
    imagesc(rot90(squeeze(max(frame2,[],3))))
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    M3 = getframe(gca);
    
%     figure(figfull)
%     fullim = [M2.cdata zeros(size(M2.cdata,1),size(M1.cdata,2)+size(M3.cdata,2)-size(M2.cdata,2),3) ; M3.cdata M1.cdata];
%     imagesc(fullim)
%     axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%     set(gca,'Position',[0 0 1 1]);
%     writeVideo(vidObj,getframe(gca));
end

%% Stitched
vidObj1 = VideoWriter('outputfull.avi','Uncompressed AVI');
vidObj1.FrameRate = m.framerate;
open(vidObj1)
fig1 = figure('Position',[0 0 ss(1)*2+ss(2)*2 (ss(3)+ss(2))*2],'Color','k');
for i = 1:ss(4)
    frame = reshape(W_final*bsxfun(@times,cmap,H_final(:,i)/max(H_final(:))),[ss(1:3) 3]);
    frame = frame*1.5;
    figure(fig1)
    imagesc([rot90(squeeze(max(frame,[],2))) zeros(ss(3),ss(2),3) ; rot90(squeeze(max(frame,[],3))) rot90(squeeze(max(frame,[],1)),2) zeros(ss(2),abs(ss(2)-ss(3)),3) ])
    axis image; axis off; set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    set(gca,'Position',[0 0 1 1]);
    writeVideo(vidObj1,getframe(fig1));
end


