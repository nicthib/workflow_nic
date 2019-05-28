% Experiment with MIDI panning settings. Looks at spatials and determines
% where they lie on the x plane, assigning a pan value to their relative
% distribution.
function panning = getpanning(W,BW,BW_u)
BW(isnan(BW)) = 0; BW_u(isnan(BW_u)) = 0;
BWl = BW_u; BWr = BW-BWl; 
for i = 1:size(W,2)
    lsum = sum(reshape(W(:,i),[sqrt(size(W,1)) sqrt(size(W,1))]).*BWl);
    rsum = sum(reshape(W(:,i),[sqrt(size(W,1)) sqrt(size(W,1))]).*BWr);
    if rsum == 0
        panning(i) = 0;
    elseif lsum == 0
        panning(i) = 127;
    else
        panning(i) = round(127*rsum/(lsum+rsum));
    end
end