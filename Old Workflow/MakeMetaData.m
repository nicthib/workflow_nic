function MakeMetaData(animalfolder)

[m.CCDpath,m.stimlist,m.pathofmetadata,m.STIMpath]=GetStimListNic(animalfolder,pathofsave);
for i=1:length(m.stimlist);
    [saved]=SaveMetaData(camera,m.stimlist{i},m.CCDpath{i},m.pathofmetadata,1,m.STIMpath{i});
end






