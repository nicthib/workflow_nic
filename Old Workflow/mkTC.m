function TCout = mkTC(data,x,y)
TCout = squeeze(squeeze(nanmean(nanmean(data(y(1):y(2),x(1):x(2),:),1),2)));