close all
piclocs = [6e6 8.31e6 1.063e7];
i = 3;
tmp = Y(piclocs(i):piclocs(i)+2.24e6,1);
%%
tmp = Y(:,1);
[PKShigh,LOCShigh]= findpeaks(diff(tmp),'MinPeakDistance',3000,'MinPeakHeight',.03);

pic = zeros(4000,numel(PKShigh));
a = 1;

for i = 1:numel(PKShigh)-1
    tmprow = tmp(LOCShigh(i):LOCShigh(i+1));
    numel(tmprow);
    if 3200 < numel(tmprow) && numel(tmprow) < 3310
             
        tmprow = [tmprow; zeros(4000-numel(tmprow),1)];
        pic(:,i) = tmprow;
    elseif 2995 < numel(tmprow) && numel(tmprow) < 3200
        tmprow = [zeros(100,1);tmprow];   
        tmprow = [tmprow; zeros(4000-numel(tmprow),1)];
        pic(:,i) = tmprow;
    end
    
end
%%
pic = pic(500:3000,1:512);
imagesc(-pic);
colormap gray
