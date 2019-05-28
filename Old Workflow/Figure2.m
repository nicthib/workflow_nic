%% Figure 2
% run - cm62_2runQ_stim3
% Components - 2,5,6,8
c = [2 5 8 10];
H_gcamp = H(c,:);
W_gcamp = W(:,c);
H_g = H_gcamp_l(c,:); W_g = W_gcamp_l(:,c);
H_h = H_chbt_l(c,:); W_h = W_chbt_l(:,c);

%% First, make an Image of components with white bkg.
close all
figure('Color','White','Position',[100 100 1000 1000]);
cmap = hsv(12);
a = 1; b = 1;
for i = 1:4
    subplot(3,4,a)
    imgtmp = reshape(W_g(:,i)*cmap(c(i),:),[m.sz m.sz 3]);
    imgtmp = imgtmp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]));
    imagesc(imgtmp)
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    subplot(3,4,a+1)
    
    imgtmp = reshape(W_h(:,i)*cmap(c(i),:),[m.sz m.sz 3]);
    imgtmp = imgtmp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3]));
    imagesc(imgtmp)
    axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    a = a+2; b = b+1;
end

%%
close all
subplot(121)
t = linspace(0,180,size(H_g,2));

plot(t,smooth(H_h(1,:))-smooth(H_h(1,:),101)+3*10^-5,'Color',cmap(c(1),:)); hold on
plot(t,smooth(H_h(2,:))-smooth(H_h(2,:),101)+2*10^-5,'Color',cmap(c(2),:))
plot(t,smooth(H_h(3,:))-smooth(H_h(3,:),101)+1*10^-5,'Color',cmap(c(3),:))
plot(t,smooth(H_h(4,:))-smooth(H_h(4,:),101)-1*10^-5,'Color',cmap(c(4),:))

set(gca,'YTick',[])
xlim([0 20]); ylim([-2 4]*1e-5)
xlabel('Time (sec)')
clc

subplot(122)
plot(t,H_g(1,:)+.3,'Color',cmap(c(1),:)); hold on
plot(t,H_g(2,:)+.2,'Color',cmap(c(2),:))
plot(t,H_g(3,:)+.1,'Color',cmap(c(3),:))
plot(t,H_g(4,:)-.1,'Color',cmap(c(4),:))

set(gca,'YTick',[])
xlim([0 20]); ylim([-.2 .4])
xlabel('Time (sec)')
clc


%% Other way to visualize...
close all
Hrs = resample(H_g',100,1)';
t = linspace(0,180,size(Hrs,2));
for i = 1:4
    tmp = Hrs(i,:); tmp(tmp<0) = 0;
    subplot(4,2,i*2-1)
    plot(t,tmp)
    xl = [0 20];
    xlim(xl)
    ylim([0 .05])
    
    
    subplot(4,2,i*2)
    [AX,H_1,H_2] = plotyy(t,round(noteroll(i,:)*127/max(max(noteroll))),t,round(tmp*127/max(tmp)));
    %y_lim = get(gca,'YLim');
    y_lim = [0 128];
    set(AX(1),'YLim',y_lim)
    set(AX(2),'YLim',y_lim)
    set(H_1(1),'LineWidth',3)
    
    set(H_1(1),'Color',cmap(c(i),:))
    
    set(AX(2),'ytick',[])
    set(AX,{'ycolor'},{'k';'k'})
    xlim(AX(1),xl)
    xlim(AX(2),xl)
end
%% Ensemble plot of MIDI notes
%subplot(212)
for i = 1:4
    plot(t,round(noteroll(i,:)*127/max(max(noteroll))),'LineWidth',3,'Color',cjet(c(i),:))
    hold on
end
%plot(t,nanmean(noteroll,1))
legend('A note','C note','E note','G note');
xlim([0 20])
ylim([0 128])
xlabel('Time (sec)')
ylabel('Note Strength')

%%
%Make a few example frames
close all

full_gcamp = zeros([m.sz m.sz 3]);
full_chbt = zeros([m.sz m.sz 3]);
act_reg1 = randi(2,[1 12])-1;
act_reg2 = randi(2,[1 12])-1;

for i = 2:12
    tmp_gcamp = reshape(W_gcamp_l(:,i)*rand*act_reg1(i)*cmap(i,:),[m.sz m.sz 3]);
    full_gcamp = full_gcamp+tmp_gcamp;
    tmp_chbt = reshape(W_chbt_l(:,i)*rand*act_reg2(i)*cmap(i,:),[m.sz m.sz 3]);
    full_chbt = full_chbt+tmp_chbt;
end
figure('Color','White','Position',[100 100 600 600]);
imagesc(full_gcamp+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3])))
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
hgexport(1,'Example image 4 gcamp.eps')

figure('Color','White','Position',[100 100 600 600]);
imagesc(full_chbt+~logical(repmat(imresize(BW,[m.sz m.sz]),[1 1 3])))
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
hgexport(2,'Example image 4 hemo.eps')

%% Ensemble plot of MIDI notes
%subplot(212)
for i = 1:4
    plot(t,round(noteroll(i,:)*127/max(max(noteroll))),'LineWidth',3,'Color',cmap(c(i),:))
    hold on
end
%plot(t,nanmean(noteroll,1))
legend('A note','C note','E note','G note');
xlim([0 20])
ylim([0 128])
xlabel('Time (sec)')
ylabel('Note Strength')



