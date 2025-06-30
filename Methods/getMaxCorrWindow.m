function windowLength = getMaxCorrWindow(corr,pvalue,offset,options)

arguments
    corr double
    pvalue double
    offset double 
    options.type = 'significant' % or 'max'
end

significantWindow = pvalue <= 0.05;
if sum(significantWindow) == 0
    significantWindow = ones(length(corr),1);
end

if strcmpi(options.type,'max')
    [~,windowLength] = max(corr(significantWindow));
    windowLength = windowLength + offset;
else
    windowLength = find(significantWindow) + offset;
end

end