function [resid,HW_gcamp_corr] = getresid(data,W,H,sz)
tmp = reshape(W*H,size(data));
resid = data-tmp;
resid(resid==0) = NaN;

% Correlate H*W with data
HW_gcamp_corr = zeros(sz,sz);
h = waitbar(0,'Making correlation map...');

for i = 1:sz
    for j = 1:sz
        tmpcorr = corrcoef(squeeze(data(i,j,:)),squeeze(tmp(i,j,:)));
        HW_gcamp_corr(i,j) = tmpcorr(1,2);
    end
    waitbar(i/sz,h)
end
close(h)