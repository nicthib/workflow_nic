% Inputs:
% mouse (string)
% run (string)
% frames to load (array)


function LoadFLIR(mouse,run,vf,fr,savepth)
mdir = findmousefolder(mouse);
webdir = fullfile(mdir,'webcam',run);
vidObj = VideoWriter(fullfile(savepth,[mouse '_' run '_' mat2str(vf(1)) '_' mat2str(numel(vf)) '.avi']));
vidObj.FrameRate = fr;
open(vidObj)
outwidth = 720;
chk_sz = 1000;
chk_st = 1:chk_sz:numel(vf);
% write chk_sz frame chunks, then append to video file
for c = 1:numel(chk_st)-1
    parfor k = 1:(chk_st(c+1)-chk_st(c))
        idx = vf(k+chk_st(c)-1);
        w0 = LoadFLIR(webdir,idx,0,outwidth);
        w1 = LoadFLIR(webdir,idx,1,outwidth);
        tmp = cat(2,w0,w1);
        tval = round(idx*10/fr)/10;
        txt = sprintf([mouse '_' run ' %.1f sec'],tval);
        a(k).cdata = uint8(mean(insertText(tmp,[10 10],txt,'FontSize',18),3));
        a(k).colormap = gray(256);
    end
    writeVideo(vidObj,a)
end
clear a
% Last chunk
parfor k = 1:(numel(vf)-chk_st(end)+1)
    idx = vf(k+chk_st(end)-1);
    w0 = LoadFLIR(webdir,idx,0,outwidth);
    w1 = LoadFLIR(webdir,idx,1,outwidth);
    tmp = cat(2,w0,w1);
    tval = round(idx*10/fr)/10;
    txt = sprintf([mouse '_' run ' %.1f sec'],tval);
    a(k).cdata = uint8(mean(insertText(tmp,[10 10],txt,'FontSize',18),3));
    a(k).colormap = gray(256);
end
disp('Almost done...')
writeVideo(vidObj,a)
clear a
close(vidObj)
disp('FLIR video written!')
