%% Get centroids and arrange
BWtmp = ones(512,512);
for i = 1:12
    Wtmp = reshape(W(:,i),[512 512]);
    s = regionprops(BWtmp, Wtmp, {'Centroid','WeightedCentroid'});
    seedpix(i,:) = s.WeightedCentroid;
end
[~,I] = sort(seedpix(:,2),'ascend');

%% Mos method
clear a song Hrs
fqs = [130 155 196 233]*2/3; fqs = [fqs fqs*2 fqs*4 fqs*8 fqs*16];
amp=10; fs=4096*4;  % sampling frequency
duration=45; values=0:1/fs:duration;
for i = 1:12
    a(i,:) = amp*sin(2*pi*fqs(I(i))*values);
    Hrs(i,:) = resample(H(i,:),length(a),length(H));
end

song = sum(a.*(Hrs),1);
song = song/max(song);
audiowrite([m.mouse '_audio.wav'],song,fs);