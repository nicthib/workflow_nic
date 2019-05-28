addpath /local_mount/space/enterprise/4/fiji-linux64-20141125/Fiji.app/scripts/
Miji

%%

for j=1:length(stims_to_use)
%for j=[16:21 157:159]
    %
    %if ~(sum(strcmp(stimlist(j),stimlist(stim5_indices))) || ~sum(strcmp(stimlist(j),stimlist(stim1_indices))))
    %    continue
    %end
    %
    tic
%cd([(STIMpath{stim5_indices(j)}) '/' stimlist{stim5_indices(j)}(find(stimlist{stim5_indices(j)}~='_')) '_webcam']);
[(STIMpath{stims_to_use(j)}) '/' stimlist{stims_to_use(j)}(find(stimlist{stims_to_use(j)}~='_')) '_webcam']
cd([(STIMpath{stims_to_use(j)}) '/' stimlist{stims_to_use(j)}(find(stimlist{stims_to_use(j)}~='_')) '_webcam']);
d=dir;
1
if length(d)<500
    d={d.name};
    cd(d{3});
end
fid=fopen('True_Frame_Rate.txt');
tmp=textscan(fid,'%s');
webcam_frameRate{j}=double(cell2mat(tmp{1}(4)));
clear webcam_movie
%toc
parfor i=0:length(dir)-4 % assumed cd'd to the right directory -- all images are total directory files minus 4 files that are dir structure or framerate saved
    webcam_movie(:,:,i+1)=rgb2gray(imread([num2str(i) '.jpg'])); %
end



tic
MIJ.run('Close All')
MIJ.run('CloseEverythingElse');
MIJ.createImage(uint8(imresize(webcam_movie,.5)));
MIJ.run('NicAnalysis3');
webcam_analyzed=(MIJ.getCurrentImage);
webcam_analyzed=webcam_analyzed;


webcam_analyzed1=webcam_analyzed;

for i=1:900
    bla=bwlabel(webcam_analyzed1(:,:,i));
    if max(max(bla))>1
        for k=1:max(max(bla))
            bla_nums(k)=sum(sum(bla==k));
        end
        region_to_keep=find(bla_nums==max(bla_nums))
        bla(bla~=region_to_keep)=0;
        webcam_analyzed1(:,:,i)=(bla>0)*255;
        %imagesc(webcam_analyzed1(:,:,i)); waitforbuttonpress;
    end
end
toc
webcam_analyzed_all(:,:,:,j)=webcam_analyzed1;
end