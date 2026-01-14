function [climit, cmap_use, clim_vec] = applyDMDColormap(ax, A, cmap0, vhold)
% applyDMDColormap
% Re-slices/flips a blue->white->red colormap so that 0 always maps to white.
%
% Inputs
%   ax    - axes handle
%   A     - response map (matrix)
%   cmap0 - options.colormap (Nx3), blue->white->red with white in the middle
%   vhold - holding potential used for selecting one-sided vs diverging scaling
%
% Outputs (optional)
%   climit   - max(abs(A(:)))
%   cmap_use - colormap applied to ax
%   clim_vec - clim limits applied to ax

    if isempty(ax) || ~isgraphics(ax, 'axes')
        error('applyDMDColormap: ax must be a valid axes handle.');
    end
    if size(cmap0,2) ~= 3
        error('applyDMDColormap: cmap0 must be Nx3.');
    end

    % Find the row closest to pure white [1 1 1]
    [~, mid] = min(sum((cmap0 - 1).^2, 2));

    climit = max(abs(A), [], 'all');
    if climit == 0 || ~isfinite(climit)
        climit = eps;
    end

    % Match your existing logic in analyzeDMDSearch.m:
    % vhold > 10  -> inhibitory (positive-only)
    % vhold < -50 -> excitatory (negative-only)
    % otherwise   -> mixed
    if vhold > 10
        % 0 = white, +max = blue
        clim_vec = [0 climit];
        cmap_use = flipud(cmap0(1:mid, :));      % white -> blue
    elseif vhold < -50
        % -max = red, 0 = white
        clim_vec = [-climit 0];
        cmap_use = flipud(cmap0(mid:end, :));    % red -> white
    else
        % -max = red, 0 = white, +max = blue
        clim_vec = [-climit climit];
        cmap_use = flipud(cmap0);                % red -> white -> blue
    end

    clim(ax, clim_vec);
    colormap(ax, cmap_use);
end
