function m = FindStimRuns(m)
m.s1 = []; m.s3 = []; m.s5 = []; m.s7 = [];
for i=1:length(m.CCDpath);
    m.ftu=matfile([m.pathofmetadata '/' m.metadataruns{i}]);
    if ~isempty(m.ftu.tstim)
        switch m.ftu.tstim
            case 1; m.s1 = [m.s1,i];
            case 3; m.s3 = [m.s3,i];
            case 5; m.s5 = [m.s5,i];
            case 7; m.s7 = [m.s7,i];
        end
    end
end

m.s1pre = m.s1(m.s1<=m.splitrun); m.s1post = m.s1(m.s1>m.splitrun);
m.s3pre = m.s3(m.s3<=m.splitrun); m.s3post = m.s3(m.s3>m.splitrun);
m.s5pre = m.s5(m.s5<=m.splitrun); m.s5post = m.s5(m.s5>m.splitrun);
m.s7pre = m.s7(m.s7<=m.splitrun); m.s7post = m.s7(m.s7>m.splitrun);
m.rsruns = 1:length(m.stimlist);
m.rsruns(m.rsruns([m.s1 m.s3 m.s5 m.s7])) = [];
m = rmfield(m,{'s1','s3','s5','s7'});
end