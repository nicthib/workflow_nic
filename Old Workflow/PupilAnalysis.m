function m = PupilAnalysis(m,webcamloc)
tic; close all; clear webcam_movie
if length(m.ttu) ~= 1
    for j=1:length(m.ttu)
        if webcamloc == 1
            cd([(m.STIMpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(find(m.stimlist{m.ttu(j)}~='_')) '_webcam']);
        elseif webcamloc == 2
            stimnum = regexp(m.stimlist{m.ttu(j)},'\d','Match');
            if numel(stimnum) > 1
                stimnum = strjoin(stimnum,'');
                cd([(m.CCDpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(1:4) '_webcam' stimnum]);
            else
                cd([(m.CCDpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(1:4) '_webcam' stimnum{1}]);
            end
        end
        parfor i=0:length(dir)-4
            if webcamloc == 1
                webcam_movie(:,:,i+1)=rgb2gray(imread([num2str(i) '.jpg']));
            elseif webcamloc == 2
                webcam_movie(:,:,i+1)=rgb2gray(imread([m.stimlist{m.ttu(j)}(1:4) num2str(i) '.jpg']));
            end
        end
        
        MIJ.run('Close All')
        MIJ.createImage(uint8(imresize(webcam_movie,.5)));
        MIJ.selectWindow('Import from Matlab');
        MIJ.run('Z Project...', 'start=50 stop=900 projection=[Min Intensity]');
        m.MIP(:,:,j) = (MIJ.getCurrentImage);
        MIJ.run('Close All')
        disp(['finished ' mat2str(j) ' out of ' mat2str(length(m.ttu)) ' runs.'])
    end
end

for j=1:length(m.ttu)
    if webcamloc == 1
        cd([(m.STIMpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(find(m.stimlist{m.ttu(j)}~='_')) '_webcam']);
    elseif webcamloc == 2
        stimnum = regexp(m.stimlist{m.ttu(j)},'\d','Match');
        if numel(stimnum) > 1
            stimnum = strjoin(stimnum,'');
            cd([(m.CCDpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(1:4) '_webcam' stimnum]);
        else
            cd([(m.CCDpath{m.ttu(j)}) '/' m.stimlist{m.ttu(j)}(1:4) '_webcam' stimnum{1}]);
        end
    end
    parfor i=0:length(dir)-4
        if webcamloc == 1
            webcam_movie(:,:,i+1)=rgb2gray(imread([num2str(i) '.jpg']));
        elseif webcamloc == 2
            webcam_movie(:,:,i+1)=rgb2gray(imread([m.stimlist{m.ttu(j)}(1:4) num2str(i) '.jpg']));
        end
    end
    
    MIJ.run('Close All')
    if j == 1
        if length(m.ttu) == 1
            MIJ.createImage(uint8(imresize(webcam_movie,.5)));
        else
            MIJ.createImage(uint8(squeeze(mean(m.MIP,3))));
        end
        MIJ.run('NicAnalysis3');
        %MIJ.createImage(uint8(imresize(webcam_movie,.5)));
        %MIJ.run('NicAnalysis4');
        newthresh = 1;
    else
        MIJ.createImage(uint8(imresize(webcam_movie,.5)));
        MIJ.run('NicAnalysis4');
    end
    
    w1 = (MIJ.getCurrentImage);
    w1std = sum(std(single(w1),[],3),3);
    w1stdx = sum(w1std,1);
    w1stdy = sum(w1std,2);
    rowstoremove = find(w1stdy==0);
    colstoremove = find(w1stdx==0);
    w1(rowstoremove,:,:) = [];
    w1(:,colstoremove,:) = [];
    
    %m.MIP(:,:,j) = min(w1,[],3);
    
    MIJ.run('Close All')
    if newthresh == 1
        qthresh(w1);
        f = warndlg('Thresholding...', 'Press ok when done thresholding.');
        waitfor(f);
        w1 = evalin('base','w1');
        m.thresh = evalin('base','thresh');
        newthresh = 0;
    else
        w1 = 255-w1;
        w1(w1 < m.thresh) = 0;
        w1(w1 >= m.thresh) = 255;
        w1 = logical(w1);
    end
    
    MIJ.createImage(uint16(w1));
    MIJ.run('NicBinarize');
    %MIJ.run('Analyze Particles...', 'size=100-10000 circularity=0.30-1.00 show=Masks display stack');
    w2=(MIJ.getCurrentImage);
    MIJ.run('Close All')
    for i=1:size(w1,3)
        bla=bwlabel(w2(:,:,i));
        if max(max(bla))>1
            for k=1:max(max(bla))
                bla_nums(k)=sum(sum(bla==k));
            end
            region_to_keep=find(bla_nums==max(bla_nums));
            region_to_keep = region_to_keep(1);
            clear bla_nums;
            bla(bla~=region_to_keep)=0;
            if sum(sum(bla)) > 3
                try
                    [X,Y] = find(bla>0); m.A(i,j) = numel(X); [m.R(i,j),m.C(i,j,:),m.Xb{i,j}]=ExactMinBoundCircle([X Y]);
                catch
                    disp('missed a convex hull')
                end
            end
        else
            if sum(sum(bla)) > 3
                try
                    [X,Y] = find(bla>0); m.A(i,j) = numel(X); [m.R(i,j),m.C(i,j,:),m.Xb{i,j}]=ExactMinBoundCircle([X Y]);
                catch
                    disp('missed a convex hull')
                end
            end
        end
    end
    disp(['finished ' mat2str(j) ' out of ' mat2str(length(m.ttu)) ' runs.'])
end
toc