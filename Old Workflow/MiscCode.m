%% MISC CODE
%%
for i = 1:numel(t1)
    stdTC(i) = nanstd(m.A1(i,:));
end; stdTC(end) = 0;
meanline = nanmean(m.A1,2); meanline(end) = 0;

figure
plot(t1,meanline,'LineWidth',3)
hold on
fill([t1';fliplr(t1)'],[meanline+stdTC/sqrt(numel(m.ttu));flipud(meanline-stdTC/sqrt(numel(m.ttu)))],[0 0 1],'linestyle','none','facealpha',.2);
cd('/local_mount/space/enterprise/4/Personal Folders/Nic/Data')
if exist([m.FileName '.mat'], 'file')
  save(m.FileName,'m','m','-append','-v7.3');
else
  save(m.FileName,'m','m','-v7.3');
end

%% Make figure of components
compsum = reshape(W,[512 512 size(W,2)]);
m.jetc = jet(size(W,2));
h = figure;
compstoplot = [1:12];
for i=1:size(W,2)
    %subplot(3,6,i)
    compcolor(:,:,:,i) = reshape(reshape(squeeze(compsum(:,:,i)),[512*512 1])*m.jetc(i,:),[512 512 3]);
end
compsumcolor = sum(compcolor(:,:,:,compstoplot),4);

imagesc(compsumcolor)
hold on
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);


%% TRASH CODE

close all
figure('Color','White','Position',[100 100 1000 1000]);
cmap = hsv(18);
cmap = cmap([2 6 11 14],:);
%imgfull = zeros([m.sz m.sz 3]);
%imgfull2 = zeros([m.sz m.sz]);
%se = strel('cube',5);
%border = imresize(BW-imerode(BW,se),[m.sz m.sz]);
a = 1;
for i = [2 3 4 16]
    subplot(2,2,a)
    imgtmp = reshape(Wl(:,i)*cmap(a,:),[m.sz1 m.sz2 3]);
    %imgtmp = imgtmp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]));
   % imgtmp = imgtmp/nanmean(nanmean(nanmean(imgtmp)));
   % imgfull = imgfull + imgtmp;
    %imgfull2 = imgfull2 + reshape(W(:,i),[m.sz m.sz]);
    imagesc(imgtmp*5)
    caxis([0 1])
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    hold on;
    a = a+1;
end

%%
close all
comps = [2 3 4 16]; a = 1;
for i = comps
    plot(H(i,:)+a*.75)
    hold on
    a = a+1;
end

%%

close all
jetc = hsv;
jetc(1,:) = 0;
j = [662 988 935 300];
for i = 1:4
    subplot(2,2,i)
    imagesc(squeeze(compvid(:,:,:,ceil(j(i)*m.fr/30))*brightness+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]))))
        axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]); colormap(jetc);

end