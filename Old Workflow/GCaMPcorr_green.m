function [gcampout, green] = GCaMPcorr_green(data_raw,m)

% gcamp = zeros(size(data_raw.g));
% h = waitbar(0,'Applying green regression to gcamp');
% for i = 1:m.sz
%     for j = 1:m.sz
%         P = polyfit(log(data_raw.g(i,j,:)),log(data_raw.b(i,j,:)),1);
%         gcamp(i,j,:) = real(exp(squeeze(log(data_raw.b(i,j,:))) - squeeze(polyval(P,log(data_raw.g(i,j,:))))));
%     end
%     waitbar(i/m.sz,h)
% end
% close(h)
% gcampout = 1-(gcamp./repmat(mean(gcamp,3),[1,1,size(data_raw.b,3)]));



gcamp = zeros(size(data_raw.g));
h = waitbar(0,'Applying green regression to gcamp');
gcamp = data_raw.b./data_raw.g;
gcampout = (gcamp./repmat(mean(gcamp,3),[1,1,size(data_raw.b,3)]))-1;

green = zeros(size(data_raw.g));
green = 1-(data_raw.g./repmat(mean(data_raw.g,3),[1,1,size(data_raw.g,3)]));