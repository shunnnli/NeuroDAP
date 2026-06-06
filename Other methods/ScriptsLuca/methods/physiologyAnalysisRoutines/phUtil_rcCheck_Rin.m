function Rin=phUtil_rcCheck_Rin(adWave, options)
arguments
    adWave
    options.mode char {mustBeMember(options.mode,['peak', 'last10%', 'last30%', 'last50%'])} = 'last30%'   
    options.baseline char {mustBeMember(options.baseline,['allpre', 'pre10ms', 'pre50ms', 'pre100ms'])} = 'pre100ms'   
    options.range double 
end

if ischar(adWave)
    if evalin('base', ['exist(''' adWave ''')' ])
        adWave=evalin('base', adWave);
    else 
        Rin=nan;
        return
    end
end

    hString=adWave.UserData.headerString;
    ai=adWave.UserData.ai;

    rc_pulse=phUtil_HeaderValue(hString, 'state.phys.internal.pulseString_RCCheck');
    current_clamp=phUtil_HeaderValue(hString,'state.phys.settings.currentClamp0');
    time_step=1000/phUtil_HeaderValue(hString,'state.phys.settings.inputRate', 1); % ms

    amplitude=phUtil_parsePulsePatternString(rc_pulse, 'amplitude');
    start_point=phUtil_parsePulsePatternString(rc_pulse, 'delay')/time_step+1;
    duration=phUtil_parsePulsePatternString(rc_pulse, 'pulseWidth')/time_step;
    end_point=start_point+duration;

    end_baseline=max(start_point-1, 1);
    switch options.baseline
        case 'allpre'
            start_baseline=1;
        case 'pre10ms'
            start_baseline=(start_point-10/time_step);
        case 'pre50ms'
            start_baseline=(start_point-50/time_step);
        case 'pre100ms'
            start_baseline=(start_point-100/time_step);
    end
    start_baseline=max(start_baseline, 1);
    baseline=mean(adWave.data(start_baseline:end_baseline));

    switch options.mode
        case 'peak'
            if amplitude<0 
                vv=min(adWave.data(start_point:end_point));
            else
                vv=min(adWave.data(start_point:end_point));
            end
        case 'last10%'
            pp=0.1;
            ss1=start_point+floor((1-pp)*duration);
            vv=mean(adWave.data(ss1:end_point));
        case 'last30%'
            pp=0.3;
            ss1=start_point+floor((1-pp)*duration);
            vv=mean(adWave.data(ss1:end_point));
        case 'last50%'
            pp=0.5;
            ss1=start_point+floor((1-pp)*duration);
            vv=mean(adWave.data(ss1:end_point));
    end
 %   [amplitude vv-baseline]
    if current_clamp
        Rin=1000*(vv-baseline)/amplitude;                
    else
        Rin=1000*amplitude/(vv-baseline);                
    end

end
    




