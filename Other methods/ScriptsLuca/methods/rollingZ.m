function fff = rollingZ(rSignal, windowSize)
    gg_m = movmean(rSignal, windowSize);
    gg_s = movstd(rSignal, windowSize);
    fff=(rSignal-gg_m)./gg_s;

    g_stdZeros = find(gg_s==0);
    if ~isempty(g_stdZeros)
        disp('WARNING: Found zeros in standard deviation. Setting Infs to 0');
        fff(g_stdZeros)=0;
    end
end

