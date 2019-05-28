% Not my code, I do't remember where this came from
function newMap = MC_nonlinearCmap(myColors, centerPoint, cLim, scalingIntensity)
    dataMax = cLim(2);
    dataMin = cLim(1);
    nColors = size(myColors,1);
    colorIdx = 1:nColors;
    colorIdx = colorIdx - (centerPoint-dataMin)*numel(colorIdx)/(dataMax-dataMin); % idx wrt center point
    colorIdx = scalingIntensity * colorIdx/max(abs(colorIdx));  % scale the range
    colorIdx = sign(colorIdx).*colorIdx.^2;
    colorIdx = colorIdx - min(colorIdx);
    colorIdx = colorIdx*nColors/max(colorIdx)+1;
    newMap = interp1(colorIdx, myColors, 1:nColors);
  end