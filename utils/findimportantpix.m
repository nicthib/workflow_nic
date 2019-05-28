% Searches along timeseries data for pixels/voxels that exceed a threshold
% a given number of times. Useful for extracting important regions from
% highly sparse datasets.
function binmask = findimportantpix(data,stdT,numT,info)
binmask = zeros(1,size(data,1));
fs = info.daq.scanRate;
parfor i = 1:size(data,1)
    % Filter
    TC = highpass(data(i,:),1/10,fs);
    stdTC = std(TC);
    if length(find(TC>stdTC*stdT)) > numT
        binmask(i) = 1;
    end
end
disp(sum(binmask(:)))