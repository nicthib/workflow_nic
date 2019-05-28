function m = PupilAnalysis2(m,webcamloc)
c = [236   501    82   301];
w1 = webcam_eye(c(3):c(4),c(1):c(2),:);
for i = 1:size(w1,3)
    w2(:,:,i) = medfilt2(w1(:,:,i),[5 5]);
end

w2(w2<2) = 0;
w2(w2>90) = 0;
w2(w2 >= 1) = 1;
w2(90 <= w2) = 1;

h = waitbar(0,'Doing a thing');
m.A = zeros(size(w1,3),1); m.R = zeros(size(w1,3),1); m.C = zeros(size(w1,3),2);
misshulls = 0;
for i = 1:size(w1,3)
    tmp = w2(:,:,i);
    tmp(isnan(tmp)) = 0;
    tmp = bwconvhull(tmp,'objects');
    bla=bwlabel(tmp);
    if sum(sum(bla)) > 0
        for k=1:max(max(bla))
            bla_nums(k)=sum(sum(bla==k));
        end
        
        region_to_keep=find(bla_nums==max(bla_nums));
        region_to_keep = region_to_keep(1);
        bla(bla~=region_to_keep)=0;
        try
            [X,Y] = find(bla>0); m.A(i) = numel(X); [m.R(i),m.C(i,:),~]=ExactMinBoundCircle([X Y]);
        catch
            misshulls = misshulls + 1;
            
        end
    end
    clear bla_nums
    waitbar(i/size(w1,3),h)
end
close(h)
disp(['missed ' mat2str(misshulls) ' convex hull(s)'])
end