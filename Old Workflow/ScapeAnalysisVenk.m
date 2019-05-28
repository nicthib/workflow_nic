%%
ss = size(SCAPE_data);
timePoints = 12:ss(4)-1;
topMIP = zeros(ss(1), ss(3), length(timePoints), 'uint16');
sideMIP = zeros(ss(1), ss(2), length(timePoints), 'uint16');
windowSize = 1;
win = [-1*windowSize:windowSize];
%%
h = waitbar(0);
for i = timePoints
    % Do a temporal smoothing of the current time point with a 3 time point
    % kernel (0.3 seconds)
    tmp1 = imgaussfilt(squeeze(mean(SCAPE_data(:, :, :, i+win), 4)), 3);
    % Do a temporal smoothing of 10 time points before (1 second before)
    tmp2 = imgaussfilt(squeeze(mean(SCAPE_data(:, :, :, i-10+win), 4)), 3);
    topMIP(:, :, i) = squeeze(max(tmp1-tmp2, [], 2));
    sideMIP(:, :, i) = squeeze(max(tmp1-tmp2, [], 3));
    waitbar(i/length(timePoints));
end
close(h);
disp('Finished')
%%
figure(1);
colormap gray
for i = 1:length(timePoints)
    subplot(211);
    imagesc(squeeze(topMIP(:, :, i))',[00 50]);
    daspect([ 2.96,1.38, 1])
        title(mat2str(i))

    subplot(212);
    imagesc(flipud(squeeze(sideMIP(:, :, i))'), [00 50]);
    daspect([ 1.14,1.38, 1])

    pause(0.01)
end

%
%%
ds = squeeze(mean(mean(reshape(SCAPE_data(:,:,:,:),[3,330,149,135,2,450]),5),1));
%% THIS IS WHAT WE USED TO CREATE TIFFS
figure
timePoints = 13:898;
temp = zeros(ss(1), ss(2), ss(3), length(timePoints), 'uint16');
topMIP = zeros(ss(1), ss(3), length(timePoints), 'uint16');
sideMIP = zeros(ss(1), ss(2), length(timePoints), 'uint16');
fi = timePoints(1)-1;

for i=timePoints
    % Do a temporal smoothing of the current time point with a 3 time point
    % kernel (0.3 seconds)
    %     temp = smooth3(mean(SCAPE_data(:, :, :,i+[-2:2]), 4)-mean(SCAPE_data(:,:,:,i+[-2:2]-10),4),'box',[3,3,3]);
    temp(:, :, :, i-fi) = smooth3(mean(SCAPE_data(:, :, :,i+[-2:2]), 4)-mean(SCAPE_data(:,:,:,i+[-2:2]-10),4),'box',[3,3,3]);
    topMIP(:, :, i-fi) = squeeze(max(squeeze(temp(:, :, :, i-fi)), [], 2));
    sideMIP(:, :, i-fi) = squeeze(max(squeeze(temp(:, :, :, i-fi)), [], 3));
    %    subplot(2,1,1)
    %     imagesc(squeeze(max(temp,[],2))');%, 1);
    %     subplot(2,1,2)
    %     imagesc(flipud(squeeze(max(temp,[],3))'));%, 1);
    %     % Do a temporal smoothing of 10 time points before (1 second before)
    %     tmp2 = imgaussfilt(squeeze(mean(green(:, :, :, i-10+win), 4)), 1);
    %     topMIP(:, :, i) = squeeze(max(tmp1-tmp2, [], 2));
    %     sideMIP(:, :, i) = squeeze(max(tmp1-tmp2, [], 3));
    %     waitbar(i/length(timePoints));
%     colormap gray;
    % pause(0.01)
    i
end
disp('Finished!!!');
%%
for i = 1:length(timePoints)
    imwrite(squeeze(topMIP(:, :, i)), img2save, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
        imwrite(squeeze(sideMIP(:, :, i)), img2save2, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
i
    
    
end
%%
for i = 1:length(timePoints)
    img2save = ['runE_dsf1_10StepDF_smoothed_' num2str(i) '.tif'];
    for j = 1:size(temp, 2)
%     imwrite(squeeze(topMIP(:, :, i)), img2save, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
        imwrite(squeeze(temp(:, j, :, i)), img2save, 'tif', 'Compression', 'none', 'WriteMode', 'Append');
    end
i
end
    
    
