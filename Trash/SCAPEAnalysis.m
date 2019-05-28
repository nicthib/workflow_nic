%% Quick load and bg subtract for dendritic data
clear
load('/local_mount/space/revault/revault1/Sam_SCAPE1p/mouse17_run3_1_forliam.mat')
bg = min(SCAPE_data(:,:,:,100:130),[],4);
ss = size(SCAPE_data);
SCAPE_data_bgs = SCAPE_data-repmat(bg,[1 1 1 ss(4)]);
textprogressbar('Smoothing...') % smooth with a [3 3 5] gaussian
for i = 1:ss(2)
    SCAPE_sm(:,i,:,:) = smooth3(squeeze(SCAPE_data_bgs(:,i,:,:)),'gaussian',[3 3 3]);
    textprogressbar(100*i/ss(2))
end
ss = size(SCAPE_sm);
textprogressbar('Smoothing done')
%% Take Jout and convert to H. currently 39 components, definitely more
disp(sprintf('Creating H...\n'))
for i = 1:numel(Jout)
    if size(Jout{i},1) ~= ss(1);
        Jout{i} = Jout{i}(251:550,:,:);
    end
    [x,y,z] = ind2sub(ss(1:3),find(Jout{i} == 1));
    cnts(i,:) = round([mean(x) mean(y) mean(z)]);
    tmp = zeros(ss(4),1); % centroid vector
    for j = 1:numel(x)
       tmp = tmp+squeeze(SCAPE_sm(x(j),y(j),z(j),:));
    end
    H(i,:) = tmp'/numel(x);
end
% Sort W and H along the longest dim.
[~,idx] = sort(cnts(:,1));
H_s = H(idx,:);
% Display component centroids
% scatter3(cnts(:,1),cnts(:,2),cnts(:,3),[],hsv(size(H,1))); 
% axis equal; xlim([250 550])

%% Perform NNLS
clear SCAPE_data SCAPE_data_bgs
SCAPE_sm = reshape(SCAPE_sm,[prod(ss(1:3)) ss(4)]); % Reshape into linear vector for NNLS

disp(sprintf('Performing NNLS...\n'))
parfor i = 1:prod(ss(1:3)) % Perform NNLS. parfor makes this less hellish
    W_s(i,:) = lsqnonneg(H_s',squeeze(double(SCAPE_sm(i,:)))');    
end
% From here on out, we use H_s and W_s to indicate they have been spatially
% sorted and cropped (to avoid the silent but deadly component mismatch)
%% Visualize components
figure('Position',[200 200 1000 700])
cmap = hsv(size(H_s,1));
for i = 44:size(H_s,1)
    subplot(221)
    imagesc(rot90(squeeze(max(reshape(W_s(:,i)*cmap(i,:),[ss(1:3) 3]),[],3)),1))
    axis image
    subplot(222)
    imagesc(rot90(squeeze(max(reshape(W_s(:,i)*cmap(i,:),[ss(1:3) 3]),[],1)),2))
    axis image
    subplot(223)
    imagesc(rot90(squeeze(max(reshape(W_s(:,i)*cmap(i,:),[ss(1:3) 3]),[],2)),1))
    axis image
    subplot(224)
    plot(H_s(i,:))
    title(['Component ' mat2str(i)])    
    colormap gray
    waitforbuttonpress    
end
%% Visualize Top and side with TCs (for beth)
%load('/local_mount/space/enterprise/4/Personal Folders/Nic/New Workflow/mouse17_components.mat');
cmap = hsv(size(H_s,1));
brightness = 1;
W_stacked1 = []; W_stacked2 = []; W_stacked3 = [];
% Here, I take each W component and stack them on top of each other to make
% a nice, colorful image of all the components.
for i = 1:size(H_s,1) 
    W_stacked1 = [W_stacked1;rot90(squeeze(max(reshape(W_s(:,i)*cmap(i,:),[ss(1:3) 3]),[],3)),1)];
    W_stacked2 = [W_stacked2;rot90(squeeze(max(reshape(W_s(:,i)*cmap(i,:),[ss(1:3) 3]),[],1)),2)];
    W_stacked3 = [W_stacked3;rot90(squeeze(max(reshape(W_s(:,i)*cmap(i,:),[ss(1:3) 3]),[],1)),3)];
    
end
%%
figure
imagesc([cat(2,W_stacked1,W_stacked2)]); axis image; set(gca,'XTick',[]); 
%imagesc([W_stacked3]); axis image; set(gca,'XTick',[]); 

set(gca,'YTick',[(1:ss(2):size(H_s,1)*ss(2))+round(ss(2)/2)]); % Mark component #s on fig
set(gca,'XTick',[]);
set(gca,'YTickLabel',1:size(H_s,1));
set(gcf,'Position',[0 0 300 800])
set(gca,'Position',[0 0 1 1]);
im1 = getframe(gcf);


% Doing the same for H. You have to crop at an odd place to make the
% components align.
figure
for i = 1:39 
    plot(H_s(i,:)-i*200,'Color',cmap(i,:))
    set(gca,'XTick',[]); 
    set(gca,'YTick',[]); 
    xlim([0 ss(4)])
    ylim([-7800 0])
    hold on
end
set(gcf,'Position',[0 0 300 800])
%set(gca,'Color','k')
set(gca,'Position',[0 0 1 1]);
im2 = getframe(gcf);

%% Create stack for exporting: prelim setup steps
%compstouse = 1:size(H,1);
%compstokill = [13 15 19 24];
%compstouse(compstokill) = [];
% Find min and max of full stack. I do this for the whole stack forst to
% get the proper dynamic range for conversion to uint16, then go back and
% compute each timepoint individually to preserve memory.
stack = W_s*H_s;
smin = min(stack(:));
smax = max(stack(:));
clear stack
%% Create stacks and save into tiffs.
cmap = hsv(numel(compstouse));
h = waitbar(0,'making stack...');
for i = 1:ss(4)
    % load individal timepoints and convert to rgb tiff matrix.
    stack = zeros([ss(1:3) 3]);
    for j = compstouse
        stack = stack+reshape(reshape(W_s(:,j)*(H_s(j,i)-H(1,i)),[prod(ss(1:3)) 1])*cmap(find(compstouse==j),:),[ss(1:3) 3]);
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

%% Post hoc analysis: residual
% Using H_s and W_s, reconstruct SCAPE_data_sm

SCAPE_fake = reshape(W_s*H_s,ss);
SCAPE_resid = abs(reshape(W_s*H_s,ss)-SCAPE_sm);

%%
figure('Position',[0 0 1000 1000])
for i = 1:ss(4)
    title(i)
    subplot(131)
    dat = squeeze(max(SCAPE_sm(:,:,:,i),[],2));
    imagesc(dat)
    axis image; colormap jet; colorbar
    caxis([0 400])
    dat2 = squeeze(max(SCAPE_fake(:,:,:,i),[],2));
    subplot(132)
    imagesc(dat2)
    axis image; colormap jet; colorbar
    caxis([0 400])
    
    subplot(133)
    imagesc(dat2-dat)
    axis image; colormap jet; colorbar
    caxis([-50 50]); colormap gray
     M(i) = getframe(gcf);
end
vid = VideoWriter('SCAPE_resid.avi','Uncompressed AVI'); 
vid.FrameRate = 10; 
open(vid); writeVideo(vid,M); close(vid)

%%

figure; for i = 1:39; 
    subplot(311); 
    imagesc(squeeze(max(reshape(W_s(:,i),[ss(1:3)]),[],2))'); 
    %axis image; 
    subplot(312); 
    imagesc(flipud(squeeze(max(reshape(W_s(:,i),[ss(1:3)]),[],3))'));
    colormap gray; 
    subplot(313)
    hold on
    plot(H_s(i,:)+i*20);
    waitforbuttonpress
end

%% NEW ANALYSIS
P = findimportantpix(reshape(SCAPE_sm,[prod(ss(1:3)),ss(4)]));
[a,b,c] = ind2sub(ss(1:3),find(P));
[IDX,cnts] = kmeans([a b c],100,'MaxIter',100,'Display','off','OnlinePhase','off','Distance','Euclidean');
% IDX has seed #
P_new = P;
P_new(find(P)) = IDX;

figure
for i = 1:100
    tmp = P_new;
    tmp(tmp ~= i) = 0;
    [a1,b1,c1] = ind2sub(ss(1:3),find(tmp));
    shp = alphaShape(a1,b1,c1,1);
    plot(shp,'EdgeColor','none','FaceColor',rand(1,3))
    hold on
    J{i} = tmp/i;
end

% J contains each seed region for getting TC
% Upscale J
for i = 1:numel(J)
    tmp = zeros([750 149 140]);
    [x,y,z] = ind2sub(size(tmp),find(J{i}));
    tmpH = zeros(ss(4),1);
    for j = 1:numel(x)
       tmpH = tmpH + squeeze(squeeze(squeeze(squeeze(mean(mean(mean(SCAPE_sm(x(j):x(j)+1,y(j):y(j)+1,z(j):z(j)+1,:),3),2),1)))));
    end
    H(i,:) = tmpH/numel(x);
    i
end
%% MUSIC

% Subtract background

H_final = H_s(compstouse,:)-repmat(H_s(1,:),[numel(compstouse) 1]);
m.fr = info.daq.scanRate;
MIDI_GCAMP(H_final,H_final,m)







