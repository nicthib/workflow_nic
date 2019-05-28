function [seedpix] = makeseedpix(W,numseeds,bilatopt,thresh)
W(W<thresh) = 0;
compsum = reshape(W,[sqrt(size(W,1)) sqrt(size(W,1)) size(W,2)]);
cmap = hsv(size(W,2)); sz = sqrt(size(W,1));
imfull = zeros(sz,sz,3);
randc = randperm(size(W,2),size(W,2));
for i=1:size(compsum,3)
    %subplot(3,6,i)
    %imagesc(compsum(:,:,i))
    %axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
    imfull = imfull+reshape(W(:,i)*cmap(randc(i),:),[sz sz 3]);
end
imagesc(imfull);
axis image; axis off; set(gca,'xcolor','r','ycolor','r','XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[]);
%%

% For ordering components (music)
if bilatopt
    q = figure;
    imagesc(sum(compsum,3))
    [~,X] = ginput(1); X = 500-round(X);
    close(q)
    seedpix = ginput(numseeds/2);
    seedpix_flip = seedpix;
    seedpix_flip(:,1) = -seedpix_flip(:,1)+750-X;
    seedpix = [seedpix;seedpix_flip];
else
    for i = 1:numseeds
        title(mat2str(i))
        seedpix(i,:) = ginput(1);
        rectangle('Position',[seedpix(i,:)-5,10,10],'EdgeColor','w')
    end
end

seedpix = round(seedpix);
[~,Order] = sort(seedpix(:,2));
seedpix = flipud(seedpix(Order,:));
close all;