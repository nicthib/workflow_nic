function BatchConvertFLIR(mouse,run,f_start,f_num,fr,pth)
path = fullfile(findmousefolder(mouse),'webcam');
webdir = fullfile(path,run);
vidObj = VideoWriter(fullfile(pth,[mouse '_' run '_' mat2str(f_start) '-' mat2str(f_start+f_num) '.avi']));
vidObj.FrameRate = fr;
open(vidObj)
vf = [1:f_num]+f_start;
outwidth = 720;
chk_sz = 1000;
chk_st = 1:chk_sz:f_num;
% write chk_sz frame chunks, then append to video file
for c = 1:numel(chk_st)-1
    parfor k = 1:(chk_st(c+1)-chk_st(c))
        w0 = LoadFLIR(webdir,vf(k+chk_st(c)-1),0,outwidth);
        w1 = LoadFLIR(webdir,vf(k+chk_st(c)-1),1,outwidth);
        tmp = cat(2,w0,w1);
        txt = sprintf([mouse '_' run ' %0.1f sec'],round(vf(k+chk_st(c)-1)/fr));
        a(k).cdata = uint8(mean(insertText(tmp,[10 10],txt,'FontSize',18),3));
        a(k).colormap = gray(256);
    end
    writeVideo(vidObj,a)
end
clear a
% Last chunk
parfor k = 1:(numel(vf)-chk_st(end)+1)
    w0 = LoadFLIR(webdir,vf(k+chk_st(end)-1),0,outwidth);
    w1 = LoadFLIR(webdir,vf(k+chk_st(end)-1),1,outwidth);
    tmp = cat(2,w0,w1);
    txt = sprintf([mouse '_' run ' %0.1f sec'],vf(k)/fr);
    a(k).cdata = uint8(mean(insertText(tmp,[10 10],txt,'FontSize',18),3));
    a(k).colormap = gray(256);
end
writeVideo(vidObj,a)
clear a
close(vidObj)
disp('FLIR video written!')
