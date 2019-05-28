% Converts continuous timecourses H (n x t) into discrete MIDI notes. 
% The m structure needs: framerate (Hz), nLEDs (integer), keys 
% (MIDI notes you want to use, 0-127), mouse run stim (optional).
function MIDImat = MIDI_GCAMP(H,m,filename)
MIDImat = [];
addpath('/local_mount/space/enterprise/4/Personal Folders/Nic/MIDI/Code')
for i = 1:size(H,1)
    Hr(i,:) = interp(H(i,:),100);
end
Hr(:,1:2) = m.thresh*1.1;
for i = 1:size(Hr,1)
    TC = Hr(i,:);
    TC(TC <= m.thresh) = NaN;
    notechange = diff(isnan(TC));
    notestart = find(notechange == -1);
    noteend = find(notechange == 1);
    if isnan(TC(1)) == 0
        notestart = [1 notestart];
    end
    if isnan(TC(end)) == 0
        noteend = [noteend numel(TC)];
    end
    
    notemag = []; notekey = [];
    for j = 1:numel(notestart)
        tmp = Hr(i,notestart(j):noteend(j));
        %ns_shift = min(find(tmp > max(tmp)/2));
        tmp(tmp < m.thresh) = 0;
        notemag(j) = max(tmp); notekey(j) = m.keys(i);
        %notestart(j) = notestart(j) + ns_shift;
    end
    MIDItmp = [notestart;noteend;notemag;notekey];
    MIDImat = [MIDImat MIDItmp];
    tmp = zeros(size(TC)); tmp(tmp == 0) = NaN;
end

MIDImat(3,:) = round(min(MIDImat(3,:)*127/max(MIDImat(3,:))+20,127));
Mfinal = zeros(size(MIDImat,2),6);
Mfinal(:,3) = MIDImat(4,:)';
Mfinal(:,4) = MIDImat(3,:);
Mfinal(:,5) = MIDImat(1,:)/(m.framerate*100/m.nLEDs)';
Mfinal(:,6) = MIDImat(2,:)/(m.framerate*100/m.nLEDs)';

midiout = matrix2midi_nic(Mfinal,300,[4,2,24,8],0);
writemidi(midiout, [filename '.mid']);
return
