function m = GetMetaData(camera,makemetadata,animalfolder,pathofsave)
[m.CCDpath,m.stimlist,m.pathofmetadata,m.STIMpath]=GetStimListNic(animalfolder,pathofsave);
slashes=strfind(m.CCDpath{1},'/');
camera = 'andor'; mouse=m.CCDpath{1}(slashes(end-2)+1:slashes(end-1)-1);
h = waitbar(0,'Making metadata...');
for i=1:length(m.stimlist);
    slashes=strfind(m.CCDpath{i},'/');
    run{i}=m.CCDpath{i}(slashes(end)+1:end);
    fid=fopen([m.CCDpath{i} '/' run{i} '_info.txt'],'r');
    txt = textscan(fid,'%s','delimiter','\t'); txt = txt{1};
    fclose(fid); clear fid;
    switch camera
        case 'andor'
            movielength(i)=str2num(txt{find(strncmpi('Movielength (s)',txt,15))+2});
            frameRate(i)=str2num(txt{find(strncmpi('frameRate (hz)',txt,14))+2});
            
            binSize(i)=str2num(txt{find(strncmpi('Bin Size',txt,8))+2});
            
            try
                LEDs{str2num(txt{strmatch('Blue',txt)+1}),i}='blue';
            end
            try
                LEDs{str2num(txt{strmatch('Green',txt)+1}),i}='green';
            end
            try
                LEDs{str2num(txt{strmatch('Red',txt)+1}),i}='red';
            end
            try
                LEDs{str2num(txt{strmatch('Cyan',txt)+1}),i}='cyan';
            end
            try
                LEDs{str2num(txt{strmatch('Speckle',txt)+1}),i}='speckle';
            end
            height(i)=str2num(txt{find(strncmpi('height',txt,6))+3});
            width(i)=str2num(txt{find(strncmpi('width',txt,5))+3});
            tpre(i)=str2num(txt{find(strncmpi('Pre-Stim (s)',txt,12))+2});
            tstim(i)=str2num(txt{find(strncmpi('Stim (s)',txt,8))+2});
            tpost(i)=str2num(txt{find(strncmpi('Post-Stim (s)',txt,13))+2});
            fname=[m.CCDpath{i} '/' m.stimlist{i} '/'];
            nFrames(i)=length(dir(fname))-2;
            if 0
                nFrames(i)=length(dir(fname))-5;
                fid1=fopen(fullfile(fname,'acquisitionmetadata.ini'));
                txt=textscan(fid1,'%s','delimiter','\t');
                txt=txt{1};
                txt=txt{end};
                nFrames(i)=nFrames(i).*str2num(txt(find(txt=='=')+1:end));
            end
            
        case 'dalsa'
            nFrames(i)=length(dir([m.CCDpath{i} '/' m.stimlist(i)]))-2;
            movielength(i)=str2num(txt{strmatch('movielength (s)',txt)+1});
            frameRate(i)=str2num(txt{strmatch('frameRate (Hz)',txt)+1});
            binSize(i)=str2num(txt{strmatch('binSize',txt)+2});
            try
                LEDs{str2num(txt{strmatch('blue',txt)+1})}(i)='blue';
            end
            try
                LEDs{str2num(txt{strmatch('green',txt)+1})}(i)='green';
            end
            try
                LEDs{str2num(txt{strmatch('red',txt)+1})}(i)='red';
            end
            try
                LEDs{str2num(txt{strmatch('cyan',txt)+1})}(i)='cyan';
            end
            try
                LEDs{str2num(txt{strmatch('speckle',txt)+1})}(i)='speckle';
            end
            height(i)=1024/binSize(i);
            width(i)=1024/binSize(i);
            tpre(i)=str2num(txt{find(strncmpi('pre_stim (s)',txt,12))+1});
            tstim(i)=str2num(txt{find(strncmpi('stim (s)',txt,8))+1});
            tpost(i)=str2num(txt{find(strncmpi('post_stim (s)',txt,13))+1});
    end
    waitbar(i/length(m.stimlist),h)
end
% get numLEDS
m.nLEDs = zeros(3,length(m.stimlist));
for i = 1:length(m.stimlist)
    for j = 1:3
        m.nLEDs(j,i) = ischar(LEDs{j,i});
    end
end
m.nLEDs = sum(m.nLEDs,1);
close(h)
% Send variables to m for output
m. Height = height; m.Width = width; m.tpre = tpre; m.tstim = tstim; m.tpost = tpost;
m.LEDs = LEDs; m.nFrames = nFrames; m.frameRate = frameRate; m.binSize = binSize; m.movielength = movielength;
