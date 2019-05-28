function [data,H,W,IDX,cntr,cntl] = preprocess_behavioral(data,m,IDX)
if m.preprocess == 1;
    disp('Preprocessing...')
    data(:,:,isnan(squeeze(nanmean(nanmean(data,2),1)))) = [];
    s = size(data); m.sz = s(1);
    disp('Rotating...')
    data = rot90(data,m.nrot);
    load(fullfile(m.savepath,[m.mouse '_mask.mat']))
    disp('Applying mask...')
    data = data.*repmat(round(imresize(BW,1/m.dsf)),[1 1 s(3)]);
    disp('Smoothing...')
    if m.highpass > 0
        data = smooth3(data,'box',[3 3 3])-smooth3(data,'box',[3 3 m.highpass]);
    else
        data = smooth3(data,'box',[3 3 3]);
    end
end
s = size(data); m.sz = s(1);

disp('Getting H from k-means...')
[H,cntl,cntr] = getHfromKmeans(data,IDX,BW,BW_u);
disp('Performing LSQ...')
[~,W] = LSQanalysis(data,H);
disp('Arranging components...')
cord = arrangecomps(H,W);
H = H(cord,:); W = W(:,cord); 
cntl = cntl(cord,:); cntr = cntr(cord,:);
IDXs = zeros(size(IDX)); a = 1;
for i = cord'
    IDXs(IDX == i) = a;
    a = a+1;
end
IDX = IDXs;
disp('DONE!')