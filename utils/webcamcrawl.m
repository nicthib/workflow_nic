function idx = webcamcrawl(path,skipn)
idx = [];
files = dir(path);
a = figure;
imnum = 3;
while imnum < numel(files)
    im = imread(fullfile(path,files(imnum).name));
    imagesc(im)
    colormap gray
    axis image
    title(strrep(files(imnum).name,'_',' '))
    keydown = waitforbuttonpress;
    if (keydown == 0)
        idx = [idx files(imnum).name];
    else
    end
    imnum = imnum+skipn;
end
close all