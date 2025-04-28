function out = preprocessSignalVCv2(in)

    Fs = 10000;
    y = lowpass(in',2000,Fs);
    yT = sgolayfilt(y,5,27); 
    out = movmedian(yT',20,2);%6

end