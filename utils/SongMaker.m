try
    H1 = m.H';
catch
    H1 = H';
end

TrialSec = 40;
FS=20000; FR = size(H1,1)/TrialSec;
for i = 1:size(H1,2)
    stdall(i) = std(H1(:,i));
end
thresh = mean(stdall);
for i = 1:size(H1,2)
   Hsc(:,i) = resample(H1(:,i),round(FS/FR),1);
   %Hsc(:,i) = Hsc(:,i)/max(Hsc(:,i));
   %Hsc(:,i) = smooth(Hsc(:,i));
end
Hsc(Hsc<0) = 0;

Ts=1/FS;
Mfinal = [];
t=1/FS:1/FS:TrialSec;

x1 = Hsc-thresh;
x2 = x1;
x2(x2 < thresh) = 0;
x2 = round(x2*127);
x2bin = sum(x2,2); x2bin(x2bin > 0) = 1;
qstart = find(diff(x2bin) == -1);
qend = find(diff(x2bin) == 1);
if qstart(1) > qend(1)
    qend(1) = [];
end
if qend(end) < qstart(end)
    qstart(end) = [];
end

qperiod = (qend-qstart)/FS;
qthresh = .5;
cctime = t(qstart(qperiod>qthresh));


%%
h = waitbar(0,'Making a MIDI!'); a = 1;
for i = 1:12%I'
    clear M idx notes tstart tend l1 r1
    
    idx = find(diff(sign(x1(:,i))));
    if x1(1,i) >= 0
        idx(1) = [];
    end
    if mod(numel(idx),2)>0
        idx(end) = [];
    end
    l1 = x1(idx,i); r1 = x1(idx+1,i);
    notes = idx + l1./(l1-r1);
    notes = [notes(1:2:end) notes(2:2:end)];
    
    N1 = size(notes,1);               % number of notes
    M1 = zeros(N1,6); 
    %nchords = size(keys,1);
    tstart = t(round(notes(:,1)));
    tend = t(round(notes(:,2))); 
    tmpM1_3 = zeros(1,N1);
   
    %for j = 1:N1
        % M1(j,3) = keys(mod(sum(find(tstart(j)>cctime)),nchords)+1,I(i));
        % ^ Alternate Method ^
        %M1(j,3) = keys(mod(sum(find(tstart(j)>cctime)),nchords)+1,a);
        %M1(j,4) = max(Hsc(round(notes(j,1)):round(notes(j,2)),i))*127+10;
    %end
    %M1(:,5) = tstart;  M1(:,6) = tend;
    
    %Mfinal = [Mfinal;M1];
    waitbar(a/i)
    a = a+1;
end


close(h)
Mfinal(:,1) = 1;
Mfinal(:,2) = 1;

%midipre = matrix2midi(Mfinal);
%writemidi(midipre, 'export.mid');


%% Making continuous MIDIS
%keys = [28 31 36 40 43 48 52 55 60 64 67 72 76 79 84 88 91 96];
% keys = [49 51 54 56 58]; keys = [keys-12 keys keys+12 keys+24];
% t = linspace(0,180,size(Hsc,1));
% h = waitbar(0,'Making you a MIDI!');
% a = 1; 
% for i = fliplr(Yrank')
%     Mfinal = [];
%     clear M1 N1 
%     notes1 = Hsc(:,i);
%     
%     N1 = length(notes1);                                  % number of notes
%     M1 = zeros(N1,6);                                       % Size of MIDI matrix
%     M1(:,3) = repmat(keys(a),[numel(N1) 1]);                % middle C = (60)
%     M1(:,5) = t; % note on
%     M1(:,6) = t+Ts;   % note off
%     M1(:,4) = Hsc(:,i)*128;
%     
%     Mfinal = [Mfinal;M1];
%     Mfinal(:,1) = 1;
%     Mfinal(:,2) = 1;
%     midi{a} = matrix2midi(Mfinal);
%     a = a+1;
%     waitbar(a/numel(Yrank))
% end
% close(h)
%%
% for i = 1:18
%     writemidi(midi{i}, ['export' mat2str(i) '.mid']);
%     i
% end


%%
% h = waitbar(0,'Making you a song!');
% clear songPost
% a = 1;
% for i = Yrank
%     t=1/FS:1/FS:TrialSec;
%     songPost(:,i) = sin(2*harmonics(a)*pi*t*basefreq);
%     songPost(:,i) = Postsc(:,i).*songPost(:,i);
%     waitbar(a/numel(Yrank))
%     a = a+1;
% end
% songPost = sum(songPost,2);
% close(h)

%%
% 
% h = waitbar(0,'Making you a song!');
% clear songPre
% a = 1;
% for i = Yrank
%     t=1/FS:1/FS:TrialSec;
%     songPre(:,i) = sin(2*harmonics(a)*pi*t*basefreq);
%     songPre(:,i) = Hsc(:,i).*songPre(:,i);
%     waitbar(a/numel(Yrank))
%     a = a+1;
% end
% songPre = sum(songPre,2);
% close(h)
