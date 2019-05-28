function webcam = LoadWebcam(WebcamPath,dsf,m,side)
% [webcam,time_std] = LoadWebcam(WebcamPath,dsf,m,side)
%
% Legacy code for PS3Eye cameras.
%
% Inputs:
% WebcamPath: full path to webcam images.
% DSF: downsample factor (Integer).
% m: Metadata from this run
% side: location of the webcam files, either 'CCD' or 'stimCCD'
%
% Outputs:
% webcam: the webcam, arranged as a 3D matrix interpolated to the framerate
% of the WFOM camera.
% time_std: a measure of the std of the webcam images over time, indicating
% movement.

textprogressbar(sprintf(['Loading ' m.run ' Webcam ']))
a = dir(fullfile(WebcamPath,'*.txt'));
FID = fopen(fullfile(WebcamPath, a.name));
tt = fgetl(FID);
frw = str2double(tt((strfind(tt,':')+1):end));
fclose(FID);
% loading
Files = dir(fullfile(WebcamPath,'*.jpg')); Files = {Files.name};
% preallocate
try
    w_raw = rgb2gray(imresize(imread(fullfile(WebcamPath,[m.run,'0.jpg'])),1/dsf));
catch
    w_raw = rgb2gray(imresize(imread(fullfile(WebcamPath,'0.jpg')),1/dsf));
end

w_raw = zeros([size(w_raw) length(Files)],'uint8');
numIm = length(Files);
try
    parfor i = 1:numIm-1;
        try
            w_raw(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[m.run num2str(i) '.jpg'])),1/dsf));
        catch
            w_raw(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[num2str(i) '.jpg'])),1/dsf));
        end
    end
    
catch
    for i = 1:numIm-1;
        try
            w_raw(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[m.run num2str(i) '.jpg'])),1/dsf));
        catch
            w_raw(:,:,i) = rgb2gray(imresize(imread(fullfile(WebcamPath,[num2str(i) '.jpg'])),1/dsf));
        end
        textprogressbar(round(i*100/numIm));
    end
end
ss = size(w_raw);
w_raw = reshape(w_raw,[ss(1)*ss(2),ss(3)]);
t1 = linspace(0,ss(3)/frw,ss(3));
t2 = linspace(0,m.nFrames/m.framerate,m.nFrames);
webcam = zeros([prod(ss(1:2)) m.nFrames],'uint8');
try
    parfor i = 1:prod(ss(1:2))
        webcam(i,:) = uint8(interp1(t1,double(w_raw(i,:)),t2));
    end
catch
    parfor i = 1:prod(ss(1:2))
        webcam(i,:) = uint8(interp1(t1,double(w_raw(i,:)),t2));
    end
end
webcam = reshape(webcam,[ss(1:2) size(webcam,2)]);


% [X,Y,Z] = meshgrid(1:ss(2),1:ss(1),linspace(0,ss(3)/frw,ss(3)));
% X = single(X); Y = single(Y); Z = single(Z);
% [X1,Y1,Z1] = meshgrid(1:ss(2),1:ss(1),linspace(0,m.nFrames/m.framerate,m.nFrames));
% X1 = single(X1); Y1 = single(Y1); Z1 = single(Z1);
% webcam = interp3(X,Y,Z,single(w_raw),X1,Y1,Z1);
%time_std = squeeze(mean(std(double(webcam(:,:,1:end-1))-double(webcam(:,:,2:end)),[],1)));
%time_std(end) = time_std(end-1);
clear tmp X* Y* Z*
%time_std = std(reshape(webcam(:,:,1:end-1)-webcam(:,:,2:end),[ss(1)*ss(2),ss(3)-1]),[],1);
textprogressbar(sprintf(' Done\n'))