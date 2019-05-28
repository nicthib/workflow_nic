function [TCout] = getwebcamrunTC(m,webcamloc)
if webcamloc == 1
    cd([(m.STIMpath{m.ttu}) '/' m.stimlist{m.ttu}(find(m.stimlist{m.ttu}~='_')) '_webcam']);
    for i=0:3:length(dir)-4
        webcam_movie(:,:,i+1)=rgb2gray(imread([num2str(i) '.jpg']));
    end
elseif webcamloc == 2
    stimnum = regexp(m.stimlist{m.ttu},'\d','Match');
    cd([(m.CCDpath{m.ttu}) '/' m.stimlist{m.ttu}(1:4) '_webcam' stimnum{1}]);
    for i=0:3:length(dir)-4
        webcam_movie(:,:,i+1)=rgb2gray(imread([m.stimlist{m.ttu}(1:4) num2str(i) '.jpg']));
    end
end
time_std_webcam=mean(std(double(webcam_movie(:,:,1:end-1))-double(webcam_movie(:,:,2:end)),[],1));
TCout=resample(squeeze((time_std_webcam)),1873,length(time_std_webcam));
mov_std=movstd(squeeze((TCout)),30);
thresh=.5+mean(TCout(find(mov_std<(3+min(mov_std)))));
