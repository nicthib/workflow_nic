% Take H and choose what you'll use as states and what you'll use as
% emissions.

StateRoll = H(1:3,:); StateRoll(isnan(StateRoll)) = 0;
EmisRoll = H(4:12,:); EmisRoll(isnan(EmisRoll)) = 0;

State = zeros(1,size(StateRoll,2));
for i = 1:size(StateRoll,1)
    State = State + StateRoll(i,:)*2^(i-1);
end

Emis = zeros(1,size(EmisRoll,2));
for i = 1:size(EmisRoll,1)
    Emis = Emis + EmisRoll(i,:)*2^(i-1);
end

%% Take the dominating component as the state, thus 12 possible states. What is a good emission? Who knows.
% normalize
for i = 1:size(H,1)
   H(i,:) = H(i,:)/max(H(i,:)); 
end
for i = 1:size(H,2)
    tmp = H(:,i);
    tmp2 = find(tmp==max(tmp));
    State(i) = tmp2(1);
end