function [wMag, wPh] = choose_sgolay_window(fHz, spanMagGHz, spanPhGHz)
    fGHz = fHz(:)/1e9;
    df   = median(diff(fGHz));
    ptsPerGHz = max(1, round(1/df));
    wMag = max(3, 2*floor((spanMagGHz*ptsPerGHz)/2)+1);
    wPh  = max(3, 2*floor((spanPhGHz *ptsPerGHz)/2)+1);
end
