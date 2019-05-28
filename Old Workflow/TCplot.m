function TCplot(m,data,opt,color)

if strcmp(opt,'ens')
    for i = 1:numel(m.ttu)
        TC_tmp = squeeze(squeeze(nanmean(nanmean(data(m.y(1,1):m.y(1,2),m.x(1,1):m.x(1,2),:,i),1),2)));
        TC_tmp(end) = NaN;
        patch(m.t,TC_tmp,'k','EdgeAlpha',.1);
        hold on
    end
elseif strcmp(opt,'std')
    stdTC = nanstd(squeeze(squeeze(nanmean(nanmean(data(m.y(1,1):m.y(1,2),m.x(1,1):m.x(1,2),:,:),1),2))),[],2);
    fill([m.t';flipud(m.t')],[mkTC(squeeze(nanmean(data,4)),m.x(1,:),...
        m.y(1,:))+stdTC;flipud(mkTC(squeeze(nanmean(data,4)),m.x(1,:),m.y(1,:))-stdTC)]...
        ,color,'linestyle','None','facealpha',.2);
    hold on
end
if nargin == 4
    plot(m.t,mkTC(squeeze(nanmean(data,4)),m.x(1,:),m.y(1,:)),'LineWidth',3,'Color',color)
else
end