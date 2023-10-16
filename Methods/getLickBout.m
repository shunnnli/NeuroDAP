function lickBout = getLickBout(lick,options)

arguments
    lick double
    options.maxILI double = 0.5 % maximum interlick interval
    options.behaviorFs double = 10000
    %options.minLicks double = 2 % min licks within a bout
    %options.minPreceedingLick double = 0.5 % min interval between preceeding lick
    %options.minProceedingLick double = 0.5 % min interval between proceeding lick
end

if isempty(lick); lickBout = []; return; end

% ILI
ILI = [5*options.behaviorFs,diff(lick)];
lickBoutIdx = find(ILI > options.maxILI*options.behaviorFs);

lickBout(:,1) = lick(lickBoutIdx)';
lickBout(:,2) = [diff(lickBoutIdx),1+(length(lick)-lickBoutIdx(end))]';


end