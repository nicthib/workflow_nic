%% Create blank midi struct
clear midib Hsc
cd(['/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/' m.mouse '/MIDI'])
%H = H_chbt_l;
for i = 1:size(H,1)
    H(i,:) = H(i,:) - smooth(H(i,:),101)';
end
H = H/max(max(H));
H(H<0) = 0;
for i = 1:size(H,1)
    Hsc(i,:) = resample(H(i,:),18000,round(m.nFrames/m.nLEDs));
    %Hsc(i,:) = smooth(Hsc(i,:))';
end
Hsc(Hsc<0) = 0;
notelength = 10; %each note takes up 10 samples,and the

% Seperate MIDI tracks need to be made because note volume is global across the
% MIDI track. j iterates through each component and writes a seperate MIDI
% file for each component.
h = waitbar(0,'Exporting MIDI files...');
compstouse = [1:12];
for j = compstouse
    clear midib
    midib = struct('format',0,'ticks_per_quarter_note',300,'track',struct);
    midib.track = struct('messages',struct);
    midib.track.messages(1).used_running_mode = [];
    midib.track.messages(1).deltatime = 0;
    midib.track.messages(1).midimeta = 0;
    midib.track.messages(1).type = 81;
    midib.track.messages(1).data = [7;161;32];
    midib.track.messages(1).chan = [];
    
    midib.track.messages(2).used_running_mode = [];
    midib.track.messages(2).deltatime = 0;
    midib.track.messages(2).midimeta = 0;
    midib.track.messages(2).type = 88;
    midib.track.messages(2).data = [4;2;24;8];
    midib.track.messages(2).chan = [];
    
    % Set note on - this remains on for the duration of the run.
    midib.track.messages(3).used_running_mode = 0;
    midib.track.messages(3).deltatime = 0;
    midib.track.messages(3).midimeta = 1;
    midib.track.messages(3).type = 144;
    midib.track.messages(3).data = [keys(j);128];
    midib.track.messages(3).chan = 1;
    
    for i = 4:numel(Hsc(1,:))+3
        midib.track.messages(i).used_running_mode = 0;
        midib.track.messages(i).deltatime = notelength*3/5;
        midib.track.messages(i).midimeta = 1;
        midib.track.messages(i).type = 176; % Message type 176 corresponds to volume modulation.
        midib.track.messages(i).data = [7;min(127,round((Hsc(j,i-3))*127))];
        midib.track.messages(i).chan = 1;
    end
    
    midib.track.messages(end+1).used_running_mode = [0];
    midib.track.messages(end).deltatime = 1000;
    midib.track.messages(end).midimeta = 1;
    midib.track.messages(end).type = 144;
    midib.track.messages(end).data = [keys(j);0];
    midib.track.messages(end).chan = 1;
    
    midib.track.messages(end+1).used_running_mode = [];
    midib.track.messages(end).deltatime = 0;
    midib.track.messages(end).midimeta = 0;
    midib.track.messages(end).type = 47;
    midib.track.messages(end).data = [];
    midib.track.messages(end).chan = [];
    
    Filename = [m.mouse '_' m.datatype];
    writemidi(midib, [Filename '_' mat2str(j) '.mid']);
    waitbar(j/numel(compstouse),h)
end
close(h)