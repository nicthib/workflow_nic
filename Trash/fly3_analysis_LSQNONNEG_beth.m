%% make low res version of SCAPE_data (2 x 2 x 2 voxels)
clear
load('/local_mount/space/revault/revault1/Sam_SCAPE1p/mouse17_run3_1_forliam.mat')
ss = size(SCAPE_data);
SCAPE_data2 = squeeze(mean(mean(mean(reshape(SCAPE_data(1:2*floor(ss(1)/2),1:2*floor(ss(2)/2),1:2*floor(ss(3)/2),:),[2,floor(ss(1)/2),2,floor(ss(2)/2),2,floor(ss(3)/2),ss(4)]),5),3),1));
%
% average over what window (to smooth, denoise)
smwin = [0:2]; %[0:2]
% subtract frame before peak
prepeak = 6; %[9]
% gaussian smoothing (spatial, xyz, must be even)
gsm = [5 5 5]; % [5 5 5]

%% Calculate magical time-course that captures each peak in the data (using the top 300 data values in the difference)
sm = squeeze(max(SCAPE_data2,[],3));
figure
i=1;
clear tempy d dd
for p = 1:4
    tempy(:,:,:,p) = smooth3(squeeze(mean(SCAPE_data(:,:,:,p+smwin),4)),'gaussian',gsm);
end
for i = 5:ss(4)
    if p==4; p = 1; end
    temp = smooth3(squeeze(mean(SCAPE_data(:,:,:,i+smwin),4)),'gaussian',gsm);
    lin = reshape(temp-tempy(:,:,:,p),[1,ss(1)*ss(2)*ss(3)]);
    p = p+1;
    [d(i,:) dd(i,:)] = hist(lin,100);
    kk=100;
    while sum(d(i,kk:end))<300
        kk=kk-1;
    end
    tcourse_sm2(i) = mean(d(i,kk:end).*dd(i,kk:end),2);
    subplot(2,1,1)
    imagesc(squeeze(max(temp-tempy(:,:,:,p),[],3)'))
    tempy(:,:,:,mod(i,4)+1) = temp;
    axis image; colormap gray
    subplot(2,1,2)
    cla
    plot(tcourse_sm2(1:i),'g')
    hold on
    pause(0.01)
end


%% Make MAPs AND Tcourses_on downsampled data for the pattern corresponding to each peak time-point (smooth data)
ss = size(SCAPE_sm2);
figure
i=1;
clear tempy
[peaks] = peakfinder(tcourse_sm2,0.1);
mapp = zeros([ss(1) ss(2) ss(3) length(peaks)]);
tcourse_samples = zeros(length(peaks),ss(4));
ii=1;
while peaks(ii)-prepeak<=0 % the first peak might be too close to the start to - prepeak
    ii=ii+1
end
m=2;
for i = peaks(ii:end)
    m=m+1;
    temp = smooth3(squeeze(mean(SCAPE_data_sm(:,:,:,i+smwin),4)),'gaussian',gsm);
    temp1 = smooth3(squeeze(mean(SCAPE_data_sm(:,:,:,i+smwin-prepeak),4)),'gaussian',gsm);
    mapp(:,:,:,m) = temp-temp1;
    subplot(4,1,1)
    imagesc(squeeze(max(mapp(:,:,:,m),[],2))');
    axis image
    subplot(4,1,2)
    imagesc(squeeze(max(mapp(:,:,:,m),[],3))');
    axis image
    colormap gray
    subplot(4,1,3)
    cla
    plot(tcourse_sm2)
    hold on
    plot(i,tcourse_sm2(i),'.r');
    
    lin = (reshape(mapp(:,:,:,m),[ss(1)*ss(2)*ss(3),1]));
    indie =  find(lin(:,1)>0.5*max(lin(:,1)));
    blank = zeros(size(lin(:,1)));
    blank(indie) = 1;
    sq = repmat(reshape(blank,[ss(1) ss(2) ss(3)]),[1,1,1,10]);
    for p = 1:10:ss(4) % this is annoying, but needed not to make a huge variable to do this masking
        tcourse_samples(m,p:p+9) = mean(mean(mean(squeeze(SCAPE_data_sm(:,:,:,p:p+9)).*sq)));
    end
    disp(i)
    subplot(4,1,4)
    plot(tcourse_samples(m,:),'-r');
    
    pause(0.1)
end
%% look at spatial mapp components (MIP)
figure
for i = 1:length(peaks)
    imagesc(squeeze(max(mapp(:,:,:,i),[],3))');
    pause(0.1)
end

%% Normalize time-courses sampled from maps
tcourse_norm = tcourse_samples./repmat(mean(tcourse_samples(:,50:100),2),[1,size(tcourse_samples,2)]);
figure; imagesc(tcourse_norm(ii:end,:)); colormap jet

%% Plot time-courses sampled together (good for seeing motion epochs)
jj = jet(length(peaks));
figure
for i = 1:length(peaks)
    plot(tcourse_norm(i,:)+0.3*i,'color',jj(i,:));
    hold on
    pause(0.1)
end

%% Correlation coefficient map of time-courses (urgh, i hate these)
for i = 1:length(peaks)
    for j = 1:length(peaks)
        x = corrcoef(tcourse_norm(i,:),tcourse_norm(j,:));
        cormap(i,j) = x(1,2);
    end
end
figure
imagesc(cormap); colormap jet

%% Eliminate periods where there is noise / motion (lazy but quick)
%% Pick regions of time-courses where motion is happening
for i = 1:10
    %or j = 1:2
    title('chose periods where motion starts and stops');
    [tttt(i,:) g] = ginput(2);
end

%% make a matrix of the parts to keep and discard

% if there is no motion:
nomo = [1:ss(4)];
mo = [];

%otherwise, run this code:
nomo = [];
mo  = []
tttt(11,1) = ss(4);
for i = 1:10
    nomo = [nomo tttt(i,2):tttt(i+1,1)];
    mo = [mo tttt(i,1):tttt(i,2)];
end


%% Looking for duplicates in peak time-courses based on their correlation with each other (didn't really use this)
figure
p=1;
used = 0;
for i = ii:length(peaks)
    %     kkk=peakfinder(corrmap(i,5:end),0.6)+4;
    clf
    kkk=find(cormap(i,:)>0.98);
    ll = length(kkk)
    ff =  ceil(ll/4);
    figure(1)
    clf
    for j = 1:ll
        subplot(ff,4,j)
        imagesc(flipud(squeeze(max(mapp(:,:,:,kkk(j)),[],3))'));
        title(kkk(j))
        colormap gray
    end
    figure(2)
    subplot(3,1,1)
    imagesc(tcourse_norm(kkk,nomo)); colormap jet
    
    if (used == i) == 0
        tnew(p,:) = mean(tcourse_norm(kkk,:),1);
        p=p+1;
        used = [used kkk];
    end
    subplot(3,1,2)
    plot(tnew(p-1,:));
    subplot(3,1,1)
    imagesc(tcourse_norm(kkk,:)); colormap jet
    pause
end


%% smooth sampled time-courses after removal of motion epochs
% smooth factor for lsqnonneg fit
smf = 5;
for i = 1:length(peaks)
    tcourse_norm_sm_mot(i,:) = smooth(tcourse_norm(i,:),smf);
end

%%
for i = 1:length(peaks)
    tnew_sm(i,:) = smooth(tnew(i,:),smf);
end

%% calculate mean intensity projection of raw data - try lsqnonneg fitting on the 2D mean before 3D (quicker)
data_flat = squeeze(mean(SCAPE_data(120:505,:,1:135,:),2));
ss = size(data_flat);
%% LSQNONNEG based on raw sampled timecourses from maps (smoothed sampled components (motion removed)). (not ICA - works) (reduced data set / mean IP)
basis = [tcourse_samples];
% amout to leave out at the front of the fit - in case of bleaching etc.
stt = 30;
cropy = 20;
neuromaps_sm = zeros([ceil(ss(1)/2)-1 ceil(ss(3)/2)-1 length(peaks)+1]);
for i = cropy:2:ss(1)-cropy-1 % selecting a sub-region of Y to speed things up
    %  for j = 1:2:ss(2)
    for k = 2:2:ss(2)-1
        there = squeeze(mean(mean(data_flat(i:i+1,k:k+1,:),2),1));
        %             for qq = 2:10 % replacing motion sections with their neighbors - this is a bad idea
        %                 there((-1+tttt(qq,1)):tttt(qq,2),1) = repmat(mean(there(tttt(qq,1)+[-5:-1],1),1),[2+tttt(qq,2)-tttt(qq,1),1]);
        %             end
        theresm = smooth(there(:),smf);
        neuromaps_sm(ceil(i/2),ceil(k/2),ii:end) = lsqnonneg(basis(ii:end,stt:end)',theresm(stt:end))';
    end
    
    disp(i)
end

%% check residuals
sn = size(neuromaps_sm);
lin = reshape(neuromaps_sm,[prod(sn(1:2)),sn(3)]);
fake = lin(:,1:end)*basis(:,1:end);
fakesq = reshape(fake,[sn(1),sn(2),ss(3)]);
mm = mean(fakesq,3);

%% compare orig mean data to Model (fake) data from lsqnonneg

%% take a look at MIPS
figure
for i = 1:length(peaks)
    subplot(2,1,1)
for i = 12:ss(3)
    
    subplot(2,1,1)
    imagesc(squeeze(data_flat(2*cropy:end-2*cropy,:,1+i)-data_flat(2*cropy:end-2*cropy,:,1+i-smf))');
    colormap gray;
    colorbar
%     c = caxis;
    caxis([0 max(c)]);
    subplot(212)
    imagesc((fakesq(cropy:end-cropy,1:end/2,i)-fakesq(cropy:end-cropy,1:end/2,i-smf))')
    colorbar
      caxis([0 max(c)]);
  pause(0.1)
    
end
    plot(tcourse_samples(i,:))
    subplot(212)
    imagesc(flipud(squeeze(max(neuromaps_sm(:,:,i),[],3))'))
    colormap gray
    title(i)
    pause(0.3)
end
%% add maps together to see spatial overlap
i=1
normmap = (1/max(max(neuromaps_sm(:,:,i))))*neuromaps_sm(:,:,i);
cumu = normmap;
jett = jet(length(peaks));
figure
for i = ii:length(peaks)
    subplot(121)
    imagesc(normmap);
    colorbar
    subplot(122)
    imagesc(cumu);
    colormap(jett)
    normmap = (1/max(max(neuromaps_sm(:,:,i))))*neuromaps_sm(:,:,i);
    if isnan(max(max(normmap)))==0
        cumu = cumu+normmap;
    end
    colorbar
    pause
end


%%
SCAPE_data_sm = squeeze(mean(mean(mean(reshape(SCAPE_data(119+[1:2*floor(ss(1)/2)],1:2*floor(ss(2)/2),1:2*floor(ss(3)/2),:),[2,floor(ss(1)/2),2,floor(ss(2)/2),2,floor(ss(3)/2),ss(4)]),5),3),1));
ss = size(SCAPE_data_sm);
%%
data_flat = squeeze(mean(SCAPE_data(120:505,:,1:135,:),2));
ss = size(data_flat);
ss = size(SCAPE_data(120:505,:,1:135,:));

%% LSQNONNEG based on raw sampled timecourses from maps (smoothed sampled components (motion removed)). (not ICA - works) (reduced data set / mean IP)
% nonmo2 = [13:93 159:261 293:389];
clear neuromaps_smz_sm
basis = [tcourse_samples];%, ones(1,size(tcourse_norm,1))'];
neuromaps_smz_sm = zeros([ceil(ss(1)) ceil(ss(2)) ceil(ss(3)) 1+length(peaks)]);

for i = 1:ss(1)-1;%120:505;%cropy/2:ss(1)-cropy/2-1
    for j = 1:ss(2)-1
        for k = 1:ss(3)-1;%1:135
            theresm = squeeze(smooth(SCAPE_data_sm(i,j,k,:),smf));%-min(smooth(SCAPE_data_sm2(i,j,k,:),smf)));
            %theresm = smooth(there,smf);
            neuromaps_smz_sm(i,j,k,:) = lsqnonneg(basis(1:end,:)',theresm(:))';
        end
    end
    disp(i)
end

%% take a look at MIPS
figure
for i = 1:length(peaks)
    subplot(2,1,1)
    imagesc(squeeze(max(neuromaps_smz_sm(:,:,:,i),[],2)))
    subplot(212)
    imagesc(flipud(squeeze(max(neuromaps_smz_sm(:,:,:,i),[],3))'))
    colormap gray
    title(i)
    pause(0.3)
end


%% check residuals
sn = size(neuromaps_smz_sm);
lin = reshape(neuromaps_smz_sm,[prod(sn(1:3)),sn(4)]);
fake = lin(:,1:end)*basis(1:end,1:end);
fakesq = reshape(fake,[sn(1) sn(2) sn(3),300]);
% lin = reshape(neuromaps_smz_nonmo,[prod(sn(1:3)),sn(4)]);
% fake = lin(:,1:end)*basis(1:end,1:end);
% fakesq2 = reshape(fake,[sn(1) sn(2) sn(3),length(peaks)]);
% mm = mean(fakesq,3);

figure
% for i = 1:length(peaks)
    subplot(2,1,1)
for i = 12:300
    
    subplot(2,1,1)
    imagesc(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))')
%     imagesc(squeeze(data_flat(2*cropy:end-2*cropy,:,1+i)-data_flat(2*cropy:end-2*cropy,:,1+i-smf))');
    colormap gray;
    colorbar
     c = caxis;
    caxis([0 max(c)]);
    subplot(212)
    imagesc(squeeze(max(fakesq(cropy:end-cropy,:,1:end,i)-fakesq(cropy:end-cropy,1:end,:,i-smf),[],2))')
    colorbar
      caxis([0 max(c)]);
  pause(0.1)
    
end
%     plot(tcourse_samples(i,:))
%     subplot(212)
%     imagesc(flipud(squeeze(max(neuromaps_sm(:,:,i),[],3))'))
%     colormap gray
%     title(i)
%     pause(0.3)
% end
%% make a bg function to correct depth attenuation (rough)
bg = squeeze(mean(mean(SCAPE_data_sm(cropy/2:end-cropy/2,:,1:end-2,:),4),3))';
bg2 = flipud(repmat(mean(bg,2),[1,161]));
bg2 = bg2-min(min(bg2));
for i = 40:50
    bg2(i,:) = bg2(1,:);
end
% bg = flipud(bg/(max(max(bg))))+0.05;
% bg(1:20,:)=1;


%% fix shift one way (diagonal shift due to oblique plane)
figure
imagesc(squeeze(mean(mean(neuromaps_smz_sm(:,:,:,10:45),4),1)));%imagesc(squeeze(max(SCAPE_data_sh(1:end,:,:),[],1))');
colormap gray
[x1 y1] = ginput(2);
shps = (max(x1)-min(x1))/(max(y1)-min(y1));

%% or just use a fixed value
shps =  1.3376;%
shps =  0.7;%

%% shift neuromaps (fix)
clear neuromaps_smzd_shift;
sn=size(neuromaps_smz_sm);
for i = ([1:sn(2)])
    m=m+1;
    bba = floor(i*shps);
%      bba = floor(i*shps);
%     dd = double(squeeze(SCAPE_data(120:505,i,1:135,:)))-minval;
    %     SCAPE_data_sh_shift(:,i,[1:ss(3)]+bba,:) = SCAPE_data_sh(:,i,:,:);%interp3([1:ss(3)],[1:ss(1)]',[1:ss(4)],squeeze(SCAPE_data_sh(:,i,:,:)),[1:ss(3)]+i*shps-bba,[1:ss(1)]',[1:ss(4)],'cubic',1);
         neuromaps_smzd_shift(:,i,[1:sn(3)]+bba,:) = interp3([1:sn(3)],[1:sn(1)]',[1:sn(4)],squeeze(neuromaps_smz_sm(:,i,1:end,:)),[1:sn(3)]-i*shps+bba,[1:sn(1)]',[1:sn(4)],'cubic',1);
%    neuromaps_smzd_shift(:,i,[1:sn(3)]-bba+ceil(shps*sn(2)),:) = interp3([1:sn(3)],[1:sn(1)]',[1:sn(4)],squeeze(neuromaps_smz_sm(:,i,1:end,:)),[1:sn(3)]+i*shps-bba,[1:sn(1)]',[1:sn(4)],'cubic',1);
%     neuromaps_smzd_shift(:,i,[1:sn(3)]+bba+ceil(shps*sn(2)),:) = interp3([1:sn(3)],[1:sn(1)]',[1:sn(4)],squeeze(neuromaps_smz_sm(:,i,1:end,:)),[1:sn(3)]+i*shps+bba,[1:sn(1)]',[1:sn(4)],'cubic',1);
          gray2(:,i,[1:sn(3)]+bba,:) = interp3([1:sn(3)],[1:sn(1)]',[1:sn(4)],squeeze(SCAPE_data_sm(:,i,1:end,1:74)-100),[1:sn(3)]-i*shps+bba,[1:sn(1)]',[1:sn(4)],'cubic',1);
   disp(i)
end



%% prepare time-courses mixed with rainbow colors
jj = jet(length(peaks)+1);
jjsh = jj(round(rand(1,length(peaks)+1)*length(peaks)+1),:);

%% Make rainbow colors (doing it this way is much faster than in the main loop)
jjsh_sq = repmat(jjsh,[1,1,ss(4)]);
ii=1
clear basiscol
for p = 1:3
    basiscol(:,:,p) = basis(ii:end,:).*squeeze(jjsh_sq(1:size(basis,1),p,:));
end

%% Prepare for multiplication below (duplicates)
% sp = size(nmaps_shift_inp);
% lin = reshape(nmaps_shift_inp,[prod(sp(1:3)),sp(4)]);
sp = size(neuromaps_smzd_shift);
lin = reshape(neuromaps_smzd_shift,[prod(sp(1:3)),sp(4)]);

% lin = lin(:,1:end-3);

%% Make nice movie that recombines colorcoded fit and data
xcal=2.2; 
ycal=1.39;
zcal = 1.13;

figure
% clf
q=0;
clear M
for i = 96:ss(4) % this just discards the parts of the data when there is motion
    q=q+1;
%     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
%     for k = ([1:sn(2)])
%         bba = floor(k*shps);
%         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
%     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol(1:end,i,p)-basiscol(1:end,i-smf,p));
        tmp =  lin*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
    end
      subplot(2,2,1)
      imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
      subplot(2,2,3)
      imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
colormap gray
    gg = caxis;
%     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
%     temp = cat(1,temp,256*ones(2,193));
%     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    temprgb = (squeeze(max(rgb(8:186,:,34:113,:)/155,[],2)));
   subplot(2,2,2)
      imagesc((permute(squeeze(temprgb),[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
       temprgb = (squeeze(max(rgb(8:186,:,34:113,:)/155,[],3)));
   subplot(2,2,4)
      imagesc(flipud(permute(temprgb,[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
%  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
%     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
%     imagesc(bigtemp)
    pause;%(0.1)
    disp(i)
%     M(q) = im2frame(uint8(256*bigtemp));
end
% currFrame = getframe (gcf);
%%
smwalk(2,:) = smooth(fly6sybrunDO(2,:),12);
smwalk(1,:) = smooth(fly6sybrunDO(1,:),12);

[ppp pp] = peakfinder(sigg)
ee=200
%% color movie with tracking etc
xcal=2.2; 
ycal=1.39;
zcal = 1.13;
ghost = mean(gray2,4)*0.03;
figure
m=1
% clf
q=0;
clear M
for i = 6:300 % this just discards the parts of the data when there is motion
    q=q+1;
    %     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
    %     for k = ([1:sn(2)])
    %         bba = floor(k*shps);
    %         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
    %     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol(1:end,i,p)-basiscol(1:end,i-smf,p));
        tmp =  lin*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
     rgb(:,:,:,p) = rgb(:,:,:,p)+ghost;
   end
    %       subplot(2,2,1)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
    %       subplot(2,2,3)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
    % colormap gray
    gg = caxis;
    %     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
    %     temp = cat(1,temp,256*ones(2,193));
    %     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    temprgb = (squeeze(max(rgb(8:186,:,34:113,:)/155,[],2)));
    subplot(2,3,2)
    imagesc([1:ss(1)]*ycal,[1:ss(3)]*xcal,flipud(permute(squeeze(temprgb),[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
    temprgb = (squeeze(max(rgb(8:186,:,34:113,:)/155,[],3)));
    axis image
    subplot(2,3,5)
    imagesc([1:ss(1)]*ycal,[1:ss(2)]*xcal,(permute(temprgb,[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
    %  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
    %     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
    %     imagesc(bigtemp)
    axis image
    title(sprintf('t = %.1f s',i/10));
    if i>10*10&i<15*10; title(sprintf('ODOR!! t = %.1f s',i/10),'color','r','FontSize',14); end
    for k = 1:5;
        subplot(2,3,3)
        imagesc(squeeze((cam(:,:,ppp(i)+k)))); colormap gray; axis image
        m=m+1;
        axis off
        subplot(2,3,6)
        cla
%         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.','color',jjr(m*2,:)); hold on
               plot(fly6sybrunDO(1,1:m*2),fly6sybrunDO(2,1:m*2),'-k'); hold on
        plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.'); hold on
        theta = atan((smwalk(2,m*2)-smwalk(2,m*2-2))./(smwalk(1,m*2)-smwalk(1,m*2-2)));
        plot(fly6sybrunDO(1,m*2)+[0 -ee*cos(theta)],fly6sybrunDO(2,m*2)+[0 -ee*sin(theta)],'-r');%,'color',jj(m*2,:)); hold on
        axis([ -4000           0       -3000        1000]);
axis([ -4000           0       -3000        1000]);
        
        %        axis off
        %         M(m) = getframe(gcf);
    end
    subplot(1,3,1);
    cla
    for w = 1:74;
        plot([1:i]/10,tcourse_norm(w,1:i)+w*0.5,'color',jjsh(w,:)); hold on
    end
    xlabel('time (s)');
    pause(0.1)
    disp(i)
    %     M(q) = im2frame(uint8(256*bigtemp));
end
% currFrame = getframe (gcf);

%%
%% color movie with tracking etc
movieName = 'simonstry8_nooghost_again3.avi';
writeObj = VideoWriter(movieName, 'Uncompressed AVI');
writeObj.FrameRate = 50;
open(writeObj);

scd=155;
scd=120;
xcal=2.2;
ycal=1.39;
zcal = 1.13;
ghost = mean(gray2,4)*0.03;
figure
set(gcf,'Position',[  616         755        1110         487]);
set(gcf,'color','k');
m=1
% clf
q=0;
clear M
for i = 6:300 % this just discards the parts of the data when there is motion
    q=q+1;
    %     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
    %     for k = ([1:sn(2)])
    %         bba = floor(k*shps);
    %         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
    %     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
             tmp = squeeze(basiscol([3:55 57:end],i,p)-basiscol([3:55 57:end],i-smf,p));
        tmp =  lin(:,[3:55 57:end])*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
%    tmp = squeeze(basiscol(1:end,i,p)-basiscol(1:end,i-smf,p));
%         tmp =  lin*tmp;
%         rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
%         %      rgb(:,:,:,p) = rgb(:,:,:,p)+ghost;
    end
    %       subplot(2,2,1)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
    %       subplot(2,2,3)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
    % colormap gray
    gg = caxis;
    %     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
    %     temp = cat(1,temp,256*ones(2,193));
    %     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],2)));
    st = size(temprgb);
    s1= subplot(2,3,1);
    imagesc([1:st(1)]*ycal,[1:st(2)]*xcal,flipud(permute(squeeze(temprgb),[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
    hold on; plot([21 70],[145 145],'w','LineWidth',2);
    xlabel('Y');
    ylabel('X');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    if i>10*10&i<15*10; title(sprintf('ODOR'),'color','r','FontSize',14); end
    
    
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],3)));
    st = size(temprgb);
    axis image
    axis off
    s2 = subplot(2,3,4);
    imagesc([1:st(1)]*ycal,[1:st(2)]*zcal,(permute(temprgb,[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
    %  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
    %     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
    %     imagesc(bigtemp)
    xlabel('Y');
    ylabel('Z');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    axis off
    axis image
    text(11,-10,sprintf('t = %.1f s',i/10),'color','w');
    for k = 1:5;
        subplot(2,3,2)
        imagesc(squeeze((cam(:,:,ppp(i)+k)))); colormap gray; axis image
        m=m+1;
        axis off
        subplot(2,3,5)
        cla
        %         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.','color',jjr(m*2,:)); hold on
        plot(fly6sybrunDO(1,1:m*2),fly6sybrunDO(2,1:m*2),'-k'); hold on
        plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.'); hold on
        theta = atan((smwalk(2,m*2)-smwalk(2,m*2-2))./(smwalk(1,m*2)-smwalk(1,m*2-2)));
        plot(fly6sybrunDO(1,m*2)+[0 -ee*cos(theta)],fly6sybrunDO(2,m*2)+[0 -ee*sin(theta)],'-r');%,'color',jj(m*2,:)); hold on
        axis([ -4000           0       -3000        1000]);
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gca,'Color','w')
        
        axis([ -4000           0       -3000        1000]);
        
        %        axis off
        %         M(m) = getframe(gcf);
        set(s1,'Position',[0.0800    0.37   0.3    0.7]);
        set(s2,'Position',[0.0800    -0.1   0.3    0.7]);
        % set(s3,'Position',[0.0800    -0.1   0.3    0.7]);
        set(s3,'Position',[0.68    0.11   0.25    0.8]);
        A = getframe(gcf);
        writeVideo(writeObj, A);
        pause(0.1)
    end
    s3=  subplot(1,3,3);
    cla
    for w = 1:74;
        plot([1:i]/10,tcourse_norm(w,1:i)+w*0.5,'color',0.1+(jjsh(w,:))*0.9); hold on
    end
    axis([0 30 0 40]);
    xlabel('time (s)');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    disp(i)
    
    %     M(q) = im2frame(uint8(256*bigtemp));
    
    
end
% currFrame = getframe (gcf);
close(writeObj)

%%
figure; for g = 1:5; for i = 1:16; subplot(4,4,i); imagesc(squeeze(max(neuromaps_smzd_shift(:,:,:,i+(g-1)*16),[],2))); colormap gray; end; pause; end
figure; for g = 1:5; for i = 1:16; subplot(4,4,i); imagesc(squeeze(max(neuromaps_smzd_shift(:,:,:,i+(g-1)*16),[],3))); title(i+(g-1)*16); colormap gray; end; pause; end

%% color movie with tracking etc GRAY
movieName = 'simonstry8_gray_fixed3_fixscale7e_sm2.avi';
writeObj = VideoWriter(movieName, 'Uncompressed AVI');
writeObj.FrameRate = 10;
open(writeObj);
profile on
scd=155;
scd=80;
xcal=2.2;
ycal=1.39;
zcal = 1.13;
ghost = mean(gray2,4)*0.03;
%figure
    s3=  subplot(1,3,3);

set(gcf,'Position',[  616         755        1110         487]);
set(gcf,'color','k');
m=1
% clf
q=0;
clear M
for i = 6:300 % this just discards the parts of the data when there is motion
    q=q+1;
    %     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
    %     for k = ([1:sn(2)])
    %         bba = floor(k*shps);
    %         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
    %     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol([3:55 57:end],i,p)-basiscol([3:55 57:end],i-smf,p));
        tmp =  lin(:,[3:55 57:end])*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
    end
    
          %        rgb(:,:,:,2) = rgb(:,:,:,2)+ghost;
%       subplot(2,2,1)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
    %       subplot(2,2,3)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
    % colormap gray
    gg = caxis;
    %     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
    %     temp = cat(1,temp,256*ones(2,193));
    %     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],2)));
    st = size(temprgb);
    s1= subplot(2,3,1);
    imagesc([1:st(1)]*ycal,[1:st(2)]*xcal,sum(flipud(permute(squeeze(temprgb),[2 1 3])),3));% = cat(2,temprgb,256*ones(193,2,3));
   % colorbar
    hold on; plot([21 70],[145 145],'w','LineWidth',2);
    xlabel('Y');
    ylabel('X');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
     caxis([0 5])
        if i>95 & i<160; caxis auto; cc = caxis; if max(cc)>5; caxis(cc); else caxis([0 5]); end; end% 
%if i>95 & i<160; caxis auto; else caxis([0 5]); end% 
    c = caxis;
    if i>10*10&i<15*10; title(sprintf('ODOR'),'color','r','FontSize',14); end
    colormap gray
    
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],3)));
    st = size(temprgb);
    axis image
    axis off
    s2 = subplot(2,3,4);
    imagesc([1:st(1)]*ycal,[1:st(2)]*zcal,sum(permute(temprgb,[2 1 3]),3));% = cat(2,temprgb,256*ones(193,2,3));
    caxis(c)
    text(140,-20,sprintf('scale max: %.1f',c(2)),'color','w')

    %colorbar
    %  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
    %     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
    %     imagesc(bigtemp)
    xlabel('Y');
    ylabel('Z');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    axis off
    axis image
    text(11,-20,sprintf('t = %.1f s',i/10),'color','w');
    for k = 1:5;
        subplot(2,3,2)
        imagesc(squeeze((cam(1:end,1:end,ppp(i)+k)))); colormap gray; axis image
        m=m+1;
        axis off
        subplot(2,3,5)
        cla
        %         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.','color',jjr(m*2,:)); hold on
        plot(fly6sybrunDO(1,1:m*2),fly6sybrunDO(2,1:m*2),'-k'); hold on
        plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.'); hold on
        theta = atan((smwalk(2,m*2)-smwalk(2,m*2-2))./(smwalk(1,m*2)-smwalk(1,m*2-2)));
        plot(fly6sybrunDO(1,m*2)+[0 -ee*cos(theta)],fly6sybrunDO(2,m*2)+[0 -ee*sin(theta)],'-r');%,'color',jj(m*2,:)); hold on
        axis([ -4000           0       -3000        1000]);
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gca,'Color','w')
        
        axis([ -4000           0       -3000        1000]);
        
        %        axis off
        %         M(m) = getframe(gcf);
        set(s1,'Position',[0.0800    0.37   0.3    0.7]);
        set(s2,'Position',[0.0800    -0.1   0.3    0.7]);
        % set(s3,'Position',[0.0800    -0.1   0.3    0.7]);
        set(s3,'Position',[0.68    0.11   0.25    0.8]);
                %pause(0.1)
    end
    s3=  subplot(1,3,3);
    cla
    for w = 1:74;
        plot([1:i]/10,tcourse_norm(w,1:i)+w*0.5,'color',0.1+(jjsh(w,:))*0.9); hold on
    end
    axis([0 30 0 40]);
    xlabel('time (s)');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    disp(i)
    A = getframe(gcf);
        writeVideo(writeObj, A);

    %     M(q) = im2frame(uint8(256*bigtemp));
    
    
end
% currFrame = getframe (gcf);
close(writeObj)

%% again - green ghost
%% color movie with tracking etc GRAY
movieName = 'simonstry8_gray_fixed3_fixscale6.avi';
writeObj = VideoWriter(movieName, 'Uncompressed AVI');
writeObj.FrameRate = 50;
open(writeObj);

scd=155;
scd=80;
xcal=2.2;
ycal=1.39;
zcal = 1.13;
ghost = mean(gray2,4)*0.03;
%figure

set(gcf,'Position',[  616         755        1110         487]);
set(gcf,'color','k');
m=1
% clf
q=0;
clear M
for i = 6:300 % this just discards the parts of the data when there is motion
    q=q+1;
    %     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
    %     for k = ([1:sn(2)])
    %         bba = floor(k*shps);
    %         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
    %     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol([3:55 57:end],i,p)-basiscol([3:55 57:end],i-smf,p));
        tmp =  lin(:,[3:55 57:end])*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
    end
    for p = 1:3
    rgb(:,:,:,p) = sum(rgb,4);
    end
                  rgb(:,:,:,2) = rgb(:,:,:,2)+ghost;
%       subplot(2,2,1)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
    %       subplot(2,2,3)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
    % colormap gray
    gg = caxis;
    %     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
    %     temp = cat(1,temp,256*ones(2,193));
    %     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],2)));
    st = size(temprgb);
    s1= subplot(2,3,1);
    imagesc([1:st(1)]*ycal,[1:st(2)]*xcal,(flipud(permute(squeeze(temprgb),[2 1 3]))));% = cat(2,temprgb,256*ones(193,2,3));
   % colorbar
    hold on; plot([21 70],[145 145],'w','LineWidth',2);
    xlabel('Y');
    ylabel('X');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    if i>95 & i<160; caxis auto; cc = caxis; if max(cc)>5; caxis(cc); else caxis([0 5]); end; end% 
    c = caxis;
    if i>10*10&i<15*10; title(sprintf('ODOR'),'color','r','FontSize',14); end
    colormap gray
    
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],3)));
    st = size(temprgb);
    axis image
    axis off
    s2 = subplot(2,3,4);
    imagesc([1:st(1)]*ycal,[1:st(2)]*zcal,(permute(temprgb,[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
    caxis(c)
    text(140,-20,sprintf('scale max: %.1f',c(2)),'color','w')

    %colorbar
    %  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
    %     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
    %     imagesc(bigtemp)
    xlabel('Y');
    ylabel('Z');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    axis off
    axis image
    text(11,-20,sprintf('t = %.1f s',i/10),'color','w');
    for k = 1:5;
        subplot(2,3,2)
        imagesc(squeeze((cam(1:2:end,1:2:end,ppp(i)+k)))); colormap gray; axis image
        m=m+1;
        axis off
        subplot(2,3,5)
        cla
        %         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.','color',jjr(m*2,:)); hold on
        plot(fly6sybrunDO(1,1:m*2),fly6sybrunDO(2,1:m*2),'-k'); hold on
        plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.'); hold on
        theta = atan((smwalk(2,m*2)-smwalk(2,m*2-2))./(smwalk(1,m*2)-smwalk(1,m*2-2)));
        plot(fly6sybrunDO(1,m*2)+[0 -ee*cos(theta)],fly6sybrunDO(2,m*2)+[0 -ee*sin(theta)],'-r');%,'color',jj(m*2,:)); hold on
        axis([ -4000           0       -3000        1000]);
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gca,'Color','w')
        
        axis([ -4000           0       -3000        1000]);
        
        %        axis off
        %         M(m) = getframe(gcf);
        set(s1,'Position',[0.0800    0.37   0.3    0.7]);
        set(s2,'Position',[0.0800    -0.1   0.3    0.7]);
        % set(s3,'Position',[0.0800    -0.1   0.3    0.7]);
        set(s3,'Position',[0.68    0.11   0.25    0.8]);
        A = getframe(gcf);
        writeVideo(writeObj, A);
        pause(0.1)
    end
    s3=  subplot(1,3,3);
    cla
    for w = 1:74;
        plot([1:i]/10,tcourse_norm(w,1:i)+w*0.5,'color',0.1+(jjsh(w,:))*0.9); hold on
    end
    axis([0 30 0 40]);
    xlabel('time (s)');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    disp(i)
    
    %     M(q) = im2frame(uint8(256*bigtemp));
    
    
end
% currFrame = getframe (gcf);
close(writeObj)


%% color movie with tracking etc GRAY BOTH
movieName = 'simonstry8_gray_both.avi';
writeObj = VideoWriter(movieName, 'Uncompressed AVI');
writeObj.FrameRate = 50;
open(writeObj);

scd=155;
scd=80;
xcal=2.2;
ycal=1.39;
zcal = 1.13;
ghost = mean(gray2,4)*0.03;
figure
set(gcf,'Position',[  616         755        1110         487]);
set(gcf,'color','k');
m=1
% clf
q=0;
clear M
for i = 6:300 % this just discards the parts of the data when there is motion
    q=q+1;
    %     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
    %     for k = ([1:sn(2)])
    %         bba = floor(k*shps);
    %         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
    %     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol([3:55 57:end],i,p)-basiscol([3:55 57:end],i-smf,p));
        tmp =  lin(:,[3:55 57:end])*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
        %      rgb(:,:,:,p) = rgb(:,:,:,p)+ghost;
    end
    %       subplot(2,2,1)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
    %       subplot(2,2,3)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
    % colormap gray
    gg = caxis;
    %     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
    %     temp = cat(1,temp,256*ones(2,193));
    %     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],2)));
    st = size(temprgb);
    s1= subplot(2,3,1)
    imagesc([1:st(1)]*ycal,[1:st(2)]*xcal,sum(flipud(permute(squeeze(temprgb),[2 1 3])),3));% = cat(2,temprgb,256*ones(193,2,3));
   % colorbar
    hold on; plot([21 70],[145 145],'w','LineWidth',2);
    xlabel('Y');
    ylabel('X');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
      c = caxis;
    if i>10*10&i<15*10; title(sprintf('ODOR'),'color','r','FontSize',14); end
    colormap gray
    
    temprgb = (squeeze(max(rgb(8:186,:,34:96,:)/scd,[],3)));
    st = size(temprgb);
    axis image
    axis off
    s2 = subplot(2,3,4)
    imagesc([1:st(1)]*ycal,[1:st(2)]*zcal,sum(permute(temprgb,[2 1 3]),3));% = cat(2,temprgb,256*ones(193,2,3));
    caxis(c)
    text(140,-20,sprintf('scale max: %.1f',c(2)),'color','w')

    %colorbar
    %  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
    %     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
    %     imagesc(bigtemp)
    xlabel('Y');
    ylabel('Z');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    axis off
    axis image
    text(11,10,sprintf('t = %.1f s',i/10),'color','w');
    for k = 1:5;
        subplot(2,3,2)
        imagesc(squeeze((cam(1:2:end,1:2:end,ppp(i)+k)))); colormap gray; axis image
        m=m+1;
        axis off
        subplot(2,3,5)
        cla
        %         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.','color',jjr(m*2,:)); hold on
        plot(fly6sybrunDO(1,1:m*2),fly6sybrunDO(2,1:m*2),'-k'); hold on
        plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.'); hold on
        theta = atan((smwalk(2,m*2)-smwalk(2,m*2-2))./(smwalk(1,m*2)-smwalk(1,m*2-2)));
        plot(fly6sybrunDO(1,m*2)+[0 -ee*cos(theta)],fly6sybrunDO(2,m*2)+[0 -ee*sin(theta)],'-r');%,'color',jj(m*2,:)); hold on
        axis([ -4000           0       -3000        1000]);
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        set(gca,'Color','w')
        
        axis([ -4000           0       -3000        1000]);
        
        %        axis off
        %         M(m) = getframe(gcf);
        set(s1,'Position',[0.0800    0.37   0.3    0.7]);
        set(s2,'Position',[0.0800    -0.1   0.3    0.7]);
        % set(s3,'Position',[0.0800    -0.1   0.3    0.7]);
        set(s3,'Position',[0.68    0.11   0.25    0.8]);
        A = getframe(gcf);
        writeVideo(writeObj, A);
        pause(0.1)
    end
    s3=  subplot(1,3,3);
    cla
    for w = 1:74;
        plot([1:i]/10,tcourse_norm(w,1:i)+w*0.5,'color',0.1+(jjsh(w,:))*0.9); hold on
    end
    axis([0 30 0 40]);
    xlabel('time (s)');
    set(gca,'XColor','w')
    set(gca,'YColor','w')
    set(gca,'Color','k')
    disp(i)
    
    %     M(q) = im2frame(uint8(256*bigtemp));
    
    
end
% currFrame = getframe (gcf);
close(writeObj)

%% color movie with tracking etc -- layers
xcal=2.2; 
ycal=1.39;
zcal = 1.13;
ghost = mean(gray2,4)*0.03;
figure
m=1
% clf
q=0;
clear M
for i = 6:300 % this just discards the parts of the data when there is motion
    q=q+1;
    %     dathere = (smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i+[-3:3]),4),'gaussian',[3 3 3])-smooth3(mean(SCAPE_data(2*cropy:end-2*cropy,:,1:end-2,i-smf+[-3:3]),4),'gaussian',[3 3 3]));%),[],3))');
    %     for k = ([1:sn(2)])
    %         bba = floor(k*shps);
    %         dathere_shift(:,k,[1:sn(3)-1]-bba+ceil(shps*sn(2))) = interp2([1:sn(3)],[1:sn(1)]',squeeze(dathere(:,k,:)),[1:sn(3)-1]+k*shps-bba,[1:sn(1)]','cubic',1);
    %     end
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol(1:end,i,p)-basiscol(1:end,i-smf,p));
        tmp =  lin(:,[3:55 57:end])*tmp([3:55 57:end]);
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
%       rgb(:,:,:,p) = rgb(:,:,:,p)+ghost;
   end
    %       subplot(2,2,1)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],2))'))
    %       subplot(2,2,3)
    %       imagesc(flipud(squeeze(max(SCAPE_data_sm(cropy:end-cropy,:,1:end,i)-SCAPE_data_sm(cropy:end-cropy,1:end,:,i-smf),[],3))'))
    % colormap gray
    gg = caxis;
    %     temp = squeeze(max(dathere_shift(1:193,:,4:end-5),[],2))'*1/100;
    %     temp = cat(1,temp,256*ones(2,193));
    %     temp = cat(1,temp,squeeze(max(dathere_shift(1:193,1:end-3,:),[],3))'*1/100);
    %     temprgb = ;
    fr = [2 3 5 6 8 9 11 12];
    for kk = 1:8
        subplot(4,3,fr(kk))
        imagesc([1:ss(1)]*ycal,[1:ss(3)]*xcal,flipud(permute(squeeze((squeeze(mean(rgb(8:186,[1:7]+(kk-1)*7,34:113,:)/150,2)))),[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
        axis image
    end
%     temprgb = (squeeze(max(rgb(8:186,:,34:113,:)/155,[],3)));
%     subplot(2,3,5)
%     imagesc([1:ss(1)]*ycal,[1:ss(2)]*xcal,(permute(temprgb,[2 1 3])));% = cat(2,temprgb,256*ones(193,2,3));
%     %  temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
%     %     bigtemp = cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
%     %     imagesc(bigtemp)
%     axis image
%     title(sprintf('t = %.1f s',i/10));
%     if i>10*10&i<15*10; title(sprintf('ODOR!! t = %.1f s',i/10),'color','r','FontSize',14); end
%     for k = 1:5;
%         subplot(2,3,3)
%         imagesc(squeeze((cam(:,:,ppp(i)+k)))); colormap gray; axis image
%         m=m+1;
%         axis off
%         subplot(2,3,6)
%         cla
% %         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.','color',jjr(m*2,:)); hold on
%                plot(fly6sybrunDO(1,1:m*2),fly6sybrunDO(2,1:m*2),'-k'); hold on
%         plot(fly6sybrunDO(1,m*2),fly6sybrunDO(2,m*2),'.'); hold on
%         theta = atan((smwalk(2,m*2)-smwalk(2,m*2-2))./(smwalk(1,m*2)-smwalk(1,m*2-2)));
%         plot(fly6sybrunDO(1,m*2)+[0 -ee*cos(theta)],fly6sybrunDO(2,m*2)+[0 -ee*sin(theta)],'-r');%,'color',jj(m*2,:)); hold on
%         axis([ -4000           0       -3000        1000]);
% axis([ -4000           0       -3000        1000]);
%         
%         %        axis off
%         %         M(m) = getframe(gcf);
%     end
    subplot(1,3,1);
    cla
    for w = 1:74;
        plot([1:i]/10,tcourse_norm(w,1:i)+w*0.5,'color',jjsh(w,:)); hold on
    end
    xlabel('time (s)');
    pause(0.1)
    disp(i)
    %     M(q) = im2frame(uint8(256*bigtemp));
end
% currFrame = getframe (gcf);
%%
vidname = 'movie.avi';
vid = VideoWriter(vidname,'Uncompressed AVI');
set(vid,'FrameRate', 10);
open(vid);
writeVideo(vid,M2)
close(vid);


%% Make tiff stacks for 3D rendering
figure
filepath = 'F:\Beth_temp\tiffs_R2D1_part14\RGB_model_all_lessbright_nogaps\';
q=0;
clear M
scc = 1+[1:sp(2)]./sp(2);
for i = 20:1500;%nomo(27:end)
    rgb = zeros(sp(1),sp(2),sp(3), 3);
    
    for p = 1:3
        tmp = squeeze(basiscol(1:end,i,p)-basiscol(1:end,i-9,p));
        tmp =  lin*tmp;
        rgb(:,:,:,p) = rgb(:,:,:,p)+reshape(tmp,[sp(1) sp(2) sp(3)]);
    end
    temprgb = (squeeze(max(rgb(:,:,1:314,:)/25,[],2)));
    temprgb = cat(2,temprgb,256*ones(271,2,3));
    temprgb = permute(cat(2,temprgb,squeeze(max(rgb(:,:,1:314,:)/35,[],3))),[2 1 3]);
    bigtemp = temprgb;% cat(2,repmat(temp,[1,1,3]),ones(433,2,3),temprgb);
    imagesc(bigtemp)
    pause(0.1)
    disp(i)
    %      M(q) = im2frame(uint8(256*bigtemp));
    imgToSave = [filepath,'col_model_', '_t', num2str(i) '.tiff'];
    clear temp2
    for j = 1:sp(2)
        temp2 = squeeze(256*rgb(:,j,:,:)/40);
        imwrite(uint8(temp2), imgToSave, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
    end
end

