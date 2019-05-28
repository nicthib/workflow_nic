% Interpolates between two timepoints of an n x 2 dataset, adding 2 new
% timepoints in between.
function Hout = interpH(H)
Hout = zeros(size(H,1),4);
for i = 1:size(H,1)
    Hout(i,:) = linspace(H(i,1),H(i,2),4);
end
Hout = Hout(:,2:3);