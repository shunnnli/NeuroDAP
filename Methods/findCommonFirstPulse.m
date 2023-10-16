function [lowFirst,highFirst,lowFirstIdx,highFirstIdx] = ...
    findCommonFirstPulse(low_upsample,highSample,lowIdx,highIdx,scalar,lowName,highName)

% Find the first common sync pulse between two systems

maxdot = 0;
for i = 1:length(highIdx)
    for j = 1:length(lowIdx)
        l = min(length(highSample)-highIdx(i),length(low_upsample)-scalar*lowIdx(j));
        dotprod = dot(highSample(highIdx(i):highIdx(i)+l),low_upsample(scalar*lowIdx(j):scalar*lowIdx(j)+l));
        disp(['dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot),...
            ', ',highName,' pulse: ',num2str(i), ', ',lowName,' pulse: ',num2str(j)]);

        if dotprod > maxdot
            x_first_idx = i;
            y_first_idx = j;
            maxdot = dotprod;

            figure;
            plot(highSample(highIdx(i):highIdx(i)+l)); hold on
            plot(2*low_upsample(scalar*lowIdx(j):scalar*lowIdx(j)+l));
            title([highName,' pulse=',num2str(i), ', ',lowName,' pulse=',num2str(j),...
               ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
        end
    end
end

highFirstIdx = x_first_idx; highFirst = highIdx(x_first_idx);
lowFirstIdx = y_first_idx; lowFirst = lowIdx(y_first_idx); 

% figure;
% plot(highSample(highIdx(highFirstIdx):highIdx(highFirstIdx)+l)); hold on
% plot(2*low_upsample(scalar*lowIdx(lowFirstIdx):scalar*lowIdx(lowFirstIdx)+l));
% title([highName,' pulse=',num2str(highFirstIdx), ', ',lowName,' pulse=',num2str(lowFirstIdx),...
%    ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
autoArrangeFigures();

end