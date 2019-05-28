function corr_flicker(data,m)
flkcor{1} = squeeze(nanmean(nanmean(data,2),1));   
flkcor{2} = squeeze(nanmean(nanmean(data(1:end/2,1:end/2,:),2),1));   
flkcor{3} = squeeze(nanmean(nanmean(data(end/2:end,1:end/2,:),2),1));   
flkcor{4} = squeeze(nanmean(nanmean(data(1:end/2,end/2:end,:),2),1));   
flkcor{5} = squeeze(nanmean(nanmean(data(end/2:end,end/2:end,:),2),1));   

for i = 1:5
    
    flkcor{i} = flkcor{i}-flkcor{i}(1);
end