bs = 10;
for i = 1:size(seedpix,1)
    H_lsq(i,:) = squeeze(squeeze(nanmean(nanmean(gcamp(seedpix(i,2)-bs:seedpix(i,2)+bs,seedpix(i,1)-bs:seedpix(i,1)+bs,:),2),1)));
end
gcamp(isnan(gcamp)) = 0;
%%
W_lsq = zeros(512,512,size(seedpix,1));
for i = 1:512
    parfor j = 1:512
        W_lsq(i,j,:) = lsqnonneg(H_lsq',squeeze(gcamp(i,j,:)));
    end
    i
end
W_lsq = reshape(W_lsq,[512*512 size(seedpix,1)]);