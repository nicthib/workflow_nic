% Removes all NaN frames from 3D datasets
function data = cleanupdata(data)
flds = fields(data);
for i = 1:numel(flds)
    if ~strcmp(flds{i},'bkg')
        data.(flds{i})(:,:,isnan(squeeze(nanmean(nanmean(data.(flds{i}),2),1)))) = [];
        data.(flds{i}) = data.(flds{i})(:,:,1:end-5);
    end
end
