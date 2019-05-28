function RunSummary(m,datads)
close all

numsub = 5;
try
    subplot(numsub,1,5)
    tcam = 1/30:1/30:round(m.nFrames/m.fr); % Needs work
    m.A1 = hampel(m.A,10);
    m.A1(m.A1 == 0) = NaN;
    m.A1(end,:) = NaN;
    m.A1(m.A1>2000) = NaN;
    stdTC = nanstd(m.A1,[],2); stdTC(isnan(stdTC)) = 0;
    ATC = nanmean(m.A1,2); ATC(isnan(ATC)) = 0;
    fill([tcam';flipud(tcam')],[ATC+stdTC;ATC-stdTC],'b','linestyle','None','facealpha',.2);
    hold on
%     for i = 1:size(m.A1,2)
%         patch(tcam,m.A1(:,i),'Black','EdgeAlpha',.1);
%         hold on
%     end
     plot(tcam,nanmean(m.A1,2),'LineWidth',3)
     ylim([min(-ATC) max(ATC)])
    
%     subplot(numsub,1,6)
%     plot(nanmean(m.pupX,2)-nanmean(m.pupX(1,:)))
%     hold on
%     plot(nanmean(m.pupY,2)-nanmean(m.pupY(1,:)))
%     legend('X','Y')
catch
    numsub = 4;
end

subplot(numsub,1,1)
TCplot(m,datads.gcamp,'ens')
TCplot(m,datads.gcamp,'std','k')
ylim([-.05 .05])

subplot(numsub,1,2)
TCplot(m,datads.chbo,'ens')
TCplot(m,datads.chbo,'std','r')
ylim([-15 15]*1e-6)

subplot(numsub,1,3)
TCplot(m,datads.chbr,'ens')
TCplot(m,datads.chbr,'std','b')
ylim([-15 15]*1e-6)

subplot(numsub,1,4)
TCplot(m,datads.chbt,'ens')
TCplot(m,datads.chbt,'std',[0 .75 0])
ylim([-15 15]*1e-6)
