% Searches along timeseries data for pixels/voxels that exceed a threshold
% a given number of times. Useful for extracting important regions from
% highly sparse datasets.
function binmask = findimportantpix(data,stdT,numT,pbleach)
binmask = zeros(1,size(data,1));
parfor i = 1:size(data,1)
    % Filter
    stdTC = std(data(i,:));
    if length(find(data(i,:) > stdTC * stdT)) > numT
        binmask(i) = 1;
    end
end
disp(sum(binmask(:)))