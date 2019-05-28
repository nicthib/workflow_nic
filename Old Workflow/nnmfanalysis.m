% Option structure components:
% method: either LSQ or NNMF.
% kern: 1x3 array that has the kernel you want to smooth the data with
% subtr_kern: 1x3 array that has the kernel you want to smooth subtract with.
% ncomps - number of components you want (for NNMF only)
% numseeds - number of seeds to pick (during NNMF)
% seedpix - used for LSQ only as the seed pixels
% bs - box size for seeds
function [W,H] = nnmfanalysis(data,ncomps,BW)
% NNMF unsupervised
tic
% if isfield(opt,'subtr_kern')
%     disp(['Smooth subtracting using kernel: ' mat2str(opt.subtr_kern)])
%     ssdata = smooth3(data,'box',opt.kern)-smooth3(data,'box',opt.subtr_kern);
%     disp('Reshaping...')
%     reshapedata = reshape(ssdata,[size(data,1)*size(data,2) size(data,3)]);
%     reshapedata_u = reshape(ssdata.*BW,[size(data,1)*size(data,2) size(data,3)]);
% else
%     disp(['Smoothing using kernel: ' mat2str(opt.kern)])
%     ssdata = smooth3(data,'box',opt.kern);
%     disp('Reshaping...')
%     reshapedata = reshape(ssdata,[size(data,1)*size(data,2) size(data,3)]);
%     reshapedata_u = reshape(ssdata.*BW,[size(data,1)*size(data,2) size(data,3)]);
% end
% clear ssdata

disp('Reshaping...')
reshapedata = reshape(data,[size(data,1)*size(data,2) size(data,3)]);
reshapedata_u = reshape(data.*repmat(BW,[1 1 size(data,3)]),[size(data,1)*size(data,2) size(data,3)]);
reshapedata(isnan(reshapedata)) = 0; reshapedata_u(isnan(reshapedata_u)) = 0;

disp('Performing NNMF...')
[~,H] = nnmf(reshapedata_u,ncomps);
W=reshape(reshapedata,[size(data,1)*size(data,2) size(data,3)])/H;
disp(['Total time: ' mat2str(toc) ' seconds'])

% To get negatives, divide W by gcamp

%% Weighted centroid
%     tmp = squeeze(compsum_u(:,:,i));
%     C=cellfun(@(n) 1:n, num2cell(size(tmp)),'uniformoutput',0);
%     [C{:}]=ndgrid(C{:});
%     C=cellfun(@(x) x(:), C,'uniformoutput',0);
%     C=[C{:}];
%     CoM(i,:)=tmp(:).'*C/sum(tmp(:));
