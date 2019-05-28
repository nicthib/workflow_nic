function gcamp=GcampGreenCorrection(m,data_raw)
parfor i=1:prod([m.height m.width])
    %tic;
    [x,y]=ind2sub([m.height m.width],i);
    p=polyfit(data_raw.g(x,y,:),data_raw.b(x,y,:),1);
    gcamp(i,:)=exp(data_raw.b(x,y,:)-(polyval(p,data_raw.g(x,y,:))));
    %tmp=csaps(1:times,(LED1(x,y,:)),10e-4,1:times);
    %time_smooth_spline(i,:)=tmp;  
    %toc
   % parfor_progress;
end
gcamp = reshape(gcamp,(size(blue_avg)));
%gcamp = gcamp-1;
