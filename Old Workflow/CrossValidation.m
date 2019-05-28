%% Cross validation
% 1. Take random run, process random minute.
% 2. Perform NNMF, starting at some small number of components.
% 3. Apply obtained W across the rest of the data, and average the
% correlations obtain ed from residual calculations.
% Repeat with higher number of components until global minima is found.

%%
clear; cd('/local_mount/space/revault/revault2/cmdata_CCD_analysis/NNMF_Summaries/cm62_2')
load('cm62_2_params.mat'); load('cm62_2_runC_stim1_all.mat', 'm');
m.sz = 128; 

m.interval = 1:1873;

%% Load in predictor dataset.
a = 0;
for i = 16:21
    disp(['Loading run ' mat2str(i) '...'])
    clear data_raw data gcampds gcamp
    m.ttu = i;
    [m,data_raw] = LoadData(m);
    data_raw.b = data_raw.b.*repmat(BW,[1 1 size(data_raw.b,3)]);
    data_raw.r = data_raw.r.*repmat(BW,[1 1 size(data_raw.b,3)]);
    data_raw.g = data_raw.g.*repmat(BW,[1 1 size(data_raw.b,3)]);
    [data.chbo,data.chbr,data.chbt] = convert_mariel(data_raw.g,data_raw.r,'g','r',m.interval,534);
    gcamp=GcampMcCorrection(data_raw.b,data.chbr,data.chbo,m.interval,m.mua1,m.mua2);
    
    gcampds = zeros(m.sz,m.sz,1873);
    for j = 0:m.sz-1
        for k = 0:m.sz-1
            gcampds(j+1,k+1,:) = nanmean(nanmean(gcamp(j*4+1:(j+1)*4,k*4+1:(k+1)*4,:),2),1);
        end
    end
    for j = 0:2
        gcamp_ALL{a*3+j+1} = gcampds(:,:,j*624+1:(j+1)*624);
    end
    a = a+1;
end

%% Check predictor against reconstructed residuals
for i = 1:18
    gcamp_ALL{i}(isnan(gcamp_ALL{i})) = 0;
end
for n = 6:20
    [W_p,H_p] = nnmfanalysis(gcamp_ALL{1},n,imresize(BW_u,[m.sz m.sz]));
    for i = 2:18
        H = reshape(gcamp_ALL{i},[m.sz*m.sz 624])'/W_p';
        
        [~,HW_corr] = getresid(gcamp_ALL{i},W_p,H',m.sz);
        corrs(n-5,i-1) = nanmean(nanmean(HW_corr));
    end
    n
end