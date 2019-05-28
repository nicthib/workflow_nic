function [data,H,W,cntr,cntl,IDX] = preprocess_behavioral(data,m,seedpix,IDX)
if m.preprocess == 1;
    disp('Preprocessing...')
    data(:,:,isnan(squeeze(nanmean(nanmean(data,2),1)))) = [];
    s = size(data); m.sz = s(1);
    disp('Smoothing...')
    if ~isempty(m.highpass)
        data = smooth3(data,'box',m.lowpass)-smooth3(data,'box',m.highpass);
    else
        data = smooth3(data,'box',m.lowpass);
    end
end
s = size(data); m.sz = s(1);

disp('Getting H...')
[H,cntl,cntr] = getHfromKmeans(data,IDX,BW,BW_u);

disp('Performing LSQ...')
[~,W] = LSQanalysis(data,H);
disp('DONE!')