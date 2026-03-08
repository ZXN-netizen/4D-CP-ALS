function [sT, sP, sR] = alignParameters(eT, eP, eR, tT, tP, tR)
    L = length(tT); sT = zeros(1,L); sP = zeros(1,L); sR = zeros(1,L); avail = 1:length(eT);
    for idx = 1:L
        if isempty(avail), break; end
        dist = sqrt(((eT(avail)-tT(idx))/90).^2 + ((eP(avail)-tP(idx))/180).^2 + ((eR(avail)-tR(idx))/7500).^2);
        [~, minIdx] = min(dist); k = avail(minIdx);
        sT(idx) = eT(k); sP(idx) = eP(k); sR(idx) = eR(k); avail(minIdx) = []; 
    end
end