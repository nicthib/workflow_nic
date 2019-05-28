% Creates a "Neural Stream", an ensemble sinusoid representation of
% timecourses in the audible space.
function NeuralStream(H,m,filename)
H(H<0) = 0;
movielength = size(H,2)*m.nLEDs/m.framerate;
fqs = [130.8 155.56 196 233]; fqs = [fqs fqs*2 fqs*4 fqs*8 fqs*16 fqs*32];
fs = 8192;  % sampling frequency
t1 = linspace(0,movielength,size(H,2));
t2 = 0:1/fs:movielength;
for i = 1:size(H,1)
    a(i,:) = sin(2*pi*fqs(i)*t2);
    Hrs(i,:) = interp1(t1,H(i,:),t2);
end
song = sum(a.*Hrs,1);
song = song/max(song);
audiowrite([m.mouse '_' m.run '_' filename '.wav'],song,fs);
disp('Done writing Neural Stream!')
