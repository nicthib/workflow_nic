% Converts continuous timecourses H (n x t) into continuously modulating
% MIDI notes. 
% The m structure needs: movielength (in seconds), keys (the MIDI keys you
% want to use, from 0-127), mouse run and stim for filename (you can just
% make your own)
function MIDI_CHBT(H,m,filename)
addpath('/local_mount/space/enterprise/4/Personal Folders/Nic/MIDI code/Code')
%% Create blank midi struct
songlength = size(H,2)*3/m.framerate;
t1 = linspace(0,songlength,size(H,2));
t2 = 1/100:1/100:songlength;
if isfield(m,'maxH')
    H = H/m.maxH;
else
    H = H/max(H(:));
end
H(H < 0) = 0;
for i = 1:size(H,1)
    Hsc(i,:) = interp1(t1,H(i,:),t2); 
end
Hsc(:,1) = 0;
Hsc(Hsc <= 0) = 0;
notelength = 10; % each volume modulation takes up 10 ticks. 1 tick is a millisecond.

% Seperate MIDI tracks need to be made because note volume is global across the
% MIDI track. j iterates through each component and writes a seperate MIDI
% file for each component.
midiout = struct('format',1,'ticks_per_quarter_note',300,'track',struct('messages',struct));
for j = 1:size(H,1)
    % Some initialization for MIDI file
    midiout.track(j) = struct('messages',struct);
    midiout.track(j).messages(1).used_running_mode = [];
    midiout.track(j).messages(1).deltatime = 0;
    midiout.track(j).messages(1).midimeta = 0;
    midiout.track(j).messages(1).type = 81;
    midiout.track(j).messages(1).data = [7;161;32];
    midiout.track(j).messages(1).chan = [];
    
    midiout.track(j).messages(2).used_running_mode = [];
    midiout.track(j).messages(2).deltatime = 0;
    midiout.track(j).messages(2).midimeta = 0;
    midiout.track(j).messages(2).type = 88;
    midiout.track(j).messages(2).data = [4;2;24;8];
    midiout.track(j).messages(2).chan = [];
    
    % Set note on - this remains on for the duration of the run.
    midiout.track(j).messages(3).used_running_mode = 0;
    midiout.track(j).messages(3).deltatime = 0;
    midiout.track(j).messages(3).midimeta = 1;
    midiout.track(j).messages(3).type = 144;
    midiout.track(j).messages(3).data = [m.keys(j);128];
    midiout.track(j).messages(3).chan = 1;
    
    for i = 1:numel(Hsc(1,:))
        midiout.track(j).messages(i+3).used_running_mode = 0;
        midiout.track(j).messages(i+3).deltatime = notelength*3/5; % I don't know why 3/5 works
        midiout.track(j).messages(i+3).midimeta = 1;
        % Message type 176 corresponds to volume modulation.
        midiout.track(j).messages(i+3).type = 176; 
        % I boost the volume by 20%
        volume = round(Hsc(j,i)*127*1.2); if isnan(volume); volume = 0; end
        midiout.track(j).messages(i+3).data = [7;min(127,volume)];
        midiout.track(j).messages(i+3).chan = 1;
    end
    
    % Finalizing end of MIDI file 
    midiout.track(j).messages(end+1).used_running_mode = 0;
    midiout.track(j).messages(end).deltatime = 1000;
    midiout.track(j).messages(end).midimeta = 1;
    midiout.track(j).messages(end).type = 144;
    midiout.track(j).messages(end).data = [m.keys(j);0];
    midiout.track(j).messages(end).chan = 1;
    
    midiout.track(j).messages(end+1).used_running_mode = [];
    midiout.track(j).messages(end).deltatime = 0;
    midiout.track(j).messages(end).midimeta = 0;
    midiout.track(j).messages(end).type = 47;
    midiout.track(j).messages(end).data = [];
    midiout.track(j).messages(end).chan = [];
    
end
writemidi(midiout, [filename '.mid']);