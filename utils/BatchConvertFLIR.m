function BatchConvertFLIR(mouse,run,f_start,f_num,pth)
path = fullfile(findmousefolder(mouse),'webcam');
webdir = fullfile(path,run);
vidObj = VideoWriter(fullfile(pth,[mouse '_' run '_' mat2str(f_start) '-' mat2str(f_start+f_num) '.avi']));
vidObj.FrameRate = 60;
open(vidObj)
vf = [1:f_num]+f_start;
chk_sz = 1000;
chk_st = 1:chk_sz:f_num;
% write chk_sz frame chunks, then append to video file
for c = 1:numel(chk_st)-1
    parfor k = 1:(chk_st(c+1)-chk_st(c))
        w0 = LoadFLIR(webdir,vf(k+chk_st(c)-1),0,2);
        w1 = LoadFLIR(webdir,vf(k+chk_st(c)-1),1,2);
        tmp = cat(2,w0,w1);
        a(k).cdata = uint8(mean(insertText(tmp,[10 500],[mouse '_' run '_' mat2str(vf(k+chk_st(c)-1))],'FontSize',18),3));
        a(k).colormap = gray(256);
    end
    writeVideo(vidObj,a)
end
% Last chunk
parfor k = 1:(numel(vf)-chk_st(end)+1)
    w0 = LoadFLIR(webdir,vf(k+chk_st(end)-1),0,2);
    w1 = LoadFLIR(webdir,vf(k+chk_st(end)-1),1,2);
    tmp = cat(2,w0,w1);
    a(k).cdata = uint8(mean(insertText(tmp,[10 500],[mouse '_' run '_' mat2str(vf(k+chk_st(end)-1))],'FontSize',18),3));
    a(k).colormap = gray(256);
end
writeVideo(vidObj,a)
clear a
close(vidObj)
print('FLIR video written!')
