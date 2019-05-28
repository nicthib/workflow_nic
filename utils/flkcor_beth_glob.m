function [tc2,facty] = flkcor_beth_glob(tc,data,ploton)
[o1 oo] = sort(tc(:,100));
for i = 1:size(tc,2)
    tmp = squeeze(tc(oo(1:end),i))';
    [Plog(i,:) V] = polyfit([0:size(tmp,1)-1],log(tmp),1);
end
P = exp(Plog);
facty = (mean(P(:,2),1)./P(:,2))';
flkcor = squeeze(nanmean(nanmean(data,2),1));
for rr = 1:size(tc,1)
    tc2(rr,:) = tc(rr,:).*facty;
    tc3(rr,:) = tc(rr,:)./flkcor';
end
if ploton
figure
subplot(131)
plot(tc')
subplot(132)
plot(tc2')
subplot(133)
plot(tc3')
end
