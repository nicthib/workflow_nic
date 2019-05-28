function [TC] = TCmaker(data,x,y)
for i = 1:size(x,1)
    TC(i,:) = squeeze(squeeze(nanmean(nanmean(data(y(i,1):y(i,2),x(i,1):x(i,2),:),2),1)));
end




