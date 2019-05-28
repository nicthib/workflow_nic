function [CCDpath,stimlist,pathofsave,STIMpath]=GetStimListNic(animalfolder,pathofsave)
if animalfolder == 0
    animalfolder=uigetdir('/local_mount/space/revault/revault2/','Pick the CCD folder');
end
if pathofsave == 0
    pathofsave=uigetdir('/local_mount/space/revault/revault2/','Pick the analysis save folder');
end
cd(animalfolder);
ccd_dir = dir;
ccd_dir = {ccd_dir.name};
subsetted=cellfun(@(x) x(strfind(x,'run')),ccd_dir,'UniformOutput',false); %run not whisker
subsetted=cellfun(@isempty,subsetted,'UniformOutput',false);
subsetted=cell2mat(subsetted);
subsetted=ones(size(subsetted))-subsetted;
subsetted=logical(subsetted);
ccd_dir=sort(ccd_dir(subsetted));
stimlist=[]; CCDpath=[]; STIMpath=[];
for i=1:length(ccd_dir)
    clear subsetted x
    ccd_dir_s=dir(ccd_dir{i});
    ccd_dir_s = {ccd_dir_s.name};
    subsetted=cellfun(@(x) x(strfind(x,'stim')),ccd_dir_s,'UniformOutput',false); %chane stim_ to stim
    subsetted=cellfun(@isempty,subsetted,'UniformOutput',false);
    subsetted=cell2mat(subsetted);
    subsetted=ones(size(subsetted))-subsetted;
    subsetted=logical(subsetted);
    ccd_dir_s=ccd_dir_s(subsetted);
    
    subsetted=cellfun(@(x) x(strfind(x,'stim0')),ccd_dir_s,'UniformOutput',false);
    subsetted=cellfun(@isempty,subsetted,'UniformOutput',false);
    subsetted=cell2mat(subsetted);
    %subsetted=ones(size(subsetted))%-subsetted;
    subsetted=logical(subsetted);
    stimlist=[stimlist,sort(ccd_dir_s(subsetted))];
    
    whatpath=what(ccd_dir{i});
    whatpath={whatpath.path};
    CCDpath=[CCDpath,repmat(whatpath,[1 length(find(subsetted))])];
end