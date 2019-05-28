close all
fig = figure('Color','Black','Position',[100 100 1000 700]);
ha = tight_subplot(2,2);
offset = 130; scale = .4;
vidObj = VideoWriter('output.avi');
vidObj.FrameRate = m.framerate/3;
open(vidObj)
for i = 1:size(gcamp,3)
    axes(ha(3));
    imagesc(imcrop(imrotate(webcam_mouse(:,:,i*3),12),[offset offset offset+(640*scale) offset+(480*scale)]));
    axis image; axis off; colormap gray; freezeColors
    axes(ha(4));
    imagesc(webcam_eye(:,:,i))
    axis image; axis off; colormap gray; freezeColors
    axes(ha(1));
    imagesc(gcamp(:,:,i))
    caxis([-.02 .075]); axis image; axis off
    colormap gray; freezeColors
    axes(ha(2));
    imagesc(chbt(:,:,i))
    caxis([-1 1]*1e-5); axis image; axis off
    colormap jet; freezeColors
    writeVideo(vidObj,getframe(fig))
end
close(vidObj)