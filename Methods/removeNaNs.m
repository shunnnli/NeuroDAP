function varargout = removeNaNs(varargin)
% removeNaNs   Remove any entries where any input is NaN
%
%   [A1r,A2r,…] = removeNaNs(A1,A2,…)
%   removes all elements in A1, A2, … at positions where *any* Ai is NaN.
%
%   All inputs must be vectors of the same length (or empty).

    nIn = nargin;
    % make everything a column vector
    for k = 1:nIn
        varargin{k} = varargin{k}(:);
    end

    % build a mask that's true wherever *all* inputs are non-NaN
    mask = true(size(varargin{1}));
    for k = 1:nIn
        mask = mask & ~isnan(varargin{k});
    end

    % apply the mask to each input
    for k = 1:nIn
        varargout{k} = varargin{k}(mask);
    end
end
