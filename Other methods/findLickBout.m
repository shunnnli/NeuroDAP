function [bout, spontaneous] = findLickBout(thresholdSamp,cur_lick,l_lick,r_lick)
% Find subsequent licks belong to the same lick bout
% Threshold: maximum inter-lick interval (in terms of sample)

if isempty(cur_lick)
    bout = []; spontaneous = []; return;
end

% if isempty(r_lick)
%     combine_lick = sort(l_lick(l_lick>cur_lick));
% else
%     combine_lick = sort([l_lick(l_lick>cur_lick),r_lick(r_lick>cur_lick)]);
% end

combine_lick = sort([l_lick(l_lick>cur_lick),r_lick(r_lick>cur_lick)]);

if isempty(combine_lick)
    bout = []; spontaneous = []; return;
else
    % disp(diff(combine_lick));
    bout_mask = logical([1 (diff(combine_lick) < thresholdSamp)]);
    bout_end_lick = find(bout_mask==0,1)-1;
    if isempty(bout_end_lick); bout_end_lick = length(bout_mask); end
    bout = combine_lick(1:bout_end_lick);
    spontaneous = combine_lick(bout_end_lick+1:end);
end

end