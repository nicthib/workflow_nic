nComps = 12; sz = 512; nFrames = 1873;
comp = zeros([m.sz m.sz round(m.nFrames/m.nLEDs) nComps],'single');
h = waitbar(0,'Making movie components...');
for i = 1:nComps
    comp(:,:,:,i) = reshape(W(:,i)*H(i,:),[sz sz nFrames]);
    waitbar(i/nComps,h);
end; close(h)

% compvid mixes the color coded components back into a 3D matrix.
cmap = hsv(nComps);
compvid = zeros([sz sz 3 nFrames],'single');
h = waitbar(0,'Making movie frames...');
for i = 1:size(comp,3)
    compvid(:,:,:,i) = reshape(reshape(squeeze(comp(:,:,i,:)),[sz^2 nComps])*cmap,[sz sz 3]);
    waitbar(i/size(comp,3),h)
end; close(h)
clear comp