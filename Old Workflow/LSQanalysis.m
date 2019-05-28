function [H,W] = LSQanalysis(data,seedpix,bs)
for i = 1:size(seedpix,1)
    H(i,:) = squeeze(squeeze(nanmean(nanmean(data(seedpix(i,2)-...
        bs:seedpix(i,2)+bs,seedpix(i,1)-bs:seedpix(i,1)+bs,:),2),1)));
end
return
data(isnan(data)) = 0;
%H(H<0) = 0;
W = zeros(size(data,1),size(data,2),size(seedpix,1));
h = waitbar(0,'Performing LSQ...');
for i = 1:size(data,1)
    for j = 1:size(data,2)
        W(i,j,:) = lsqnonneg(H',squeeze(data(i,j,:)));
    end
    waitbar(i/size(data,1),h);
end
close(h)
W = reshape(W,[size(data,1)^2 size(seedpix,1)]);