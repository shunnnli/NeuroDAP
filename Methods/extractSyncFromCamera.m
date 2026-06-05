function [sync, intensity] = extractSyncFromCamera(aviPath, options)
    % extractSyncFromCamera Extract sync ROI intensity from a camera AVI.
    %
    %   intensity = extractSyncFromCamera(aviPath) reads the AVI file at aviPath
    %   and returns the mean grayscale intensity within the predetermined ROI for
    %   each frame.
    %
    %   intensity = extractSyncFromCamera(aviPath, userInspect=true) opens a GUI
    %   on the first frame so the ROI can be adjusted before extraction.
    %
    %   [intensity, roi] = extractSyncFromCamera(...) also returns the ROI used,
    %   formatted as [x y width height].
    %
    %   Optional name-value inputs:
    %       userInspect - logical, whether to inspect/adjust ROI interactively
    %       roi         - predetermined ROI, [x y width height]
    
    arguments
        aviPath (1,1) string
        options.inspectROI (1,1) logical = false
        options.roi (1,4) double = [310 459 50 50]
        options.saveROI logical = true
    end
    
    if ~isfile(aviPath)
        error('extractSyncFromCamera:FileNotFound', ...
            'Could not find AVI file: %s', char(aviPath));
    end
    
    % Convert to mp4 so matlab can read
    aviPath = char(strtrim(string(aviPath)));
    aviPath = osPathSwitch(aviPath);
    [aviFolder, aviName] = fileparts(aviPath);

    mp4Folder = fullfile(aviFolder, 'camera-processed');
    if ~exist(mp4Folder, 'dir')
        mkdir(mp4Folder);
    end
    mp4Path = fullfile(aviFolder,'camera-processed', [aviName '.mp4']);
    
    if ~isfile(mp4Path)
        if ispc
            ffmpegPath = 'ffmpeg';
        else
            ffmpegPath = '/opt/homebrew/bin/ffmpeg';
        end
    
        cmd = sprintf(['%s -y -i %s ' ...
            '-c:v libx264 -pix_fmt yuv420p -crf 18 -preset fast %s'], ...
            shellQuote(ffmpegPath), shellQuote(aviPath), shellQuote(mp4Path));
        
        if ispc
            [status, msg] = dos(cmd);
        else
            [status, msg] = system(cmd);
        end
        
        if status ~= 0
            error("Could not convert AVI for MATLAB VideoReader:\n%s", msg);
        end
    end
    
    % read video
    v = VideoReader(mp4Path);
    firstFrame = readFrame(v);
    frameSize = [size(firstFrame,1), size(firstFrame,2)];
    
    roi = clampRoi(options.roi, frameSize);
    
    if options.inspectROI
        roi = inspectRoi(firstFrame, roi);
        roi = clampRoi(roi, frameSize);
    end
    
    % Save ROI if necessary
    if options.saveROI
        [aviFolder, ~] = fileparts(char(aviPath));
        roiPngPath = fullfile(aviFolder, 'Camera-SyncROI.png');
        
        fig = figure('Visible', 'off', 'Color', 'w');
        imshow(firstFrame); hold on;
        rectangle('Position', roi, 'EdgeColor', 'r', 'LineWidth', 2); hold on;
        
        roiText = sprintf('ROI: [x = %.1f, y = %.1f, w = %.1f, h = %.1f]', ...
            roi(1), roi(2), roi(3), roi(4));
        
        text(0.02, 0.98, roiText,'Units','normalized','Color','w','BackgroundColor','k', ...
            'FontSize', 12, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

        exportgraphics(gca, roiPngPath, 'Resolution', 200);
        close(fig);
    end
    
    v.CurrentTime = 0;
    maxFrames = max(1, ceil(v.Duration * v.FrameRate));
    intensity = nan(maxFrames, 1);
    
    frameIdx = 0;
    disp('Ongoing: extracting sync pulse from video...');
    while hasFrame(v)
        frame = readFrame(v);
        frameIdx = frameIdx + 1;
    
        if frameIdx > numel(intensity)
            intensity = [intensity; nan(maxFrames, 1)]; %#ok<AGROW>
        end
    
        grayFrame = frameToGray(frame);
        x = roi(1);
        y = roi(2);
        w = roi(3);
        h = roi(4);
        roiPixels = grayFrame(y:y+h-1, x:x+w-1);
    
        intensity(frameIdx) = mean(roiPixels(:));
    end
    
    intensity = intensity(1:frameIdx);
    
    lowVal = prctile(intensity, 10);
    highVal = prctile(intensity, 90);
    threshold = (lowVal + highVal) / 2;
    sync = double(intensity > threshold);
    disp('Finished: extracting sync pulse from video');

end

function roi = inspectRoi(frame, roi)
    figureHandle = figure('Name', 'Inspect Sync ROI', 'NumberTitle', 'off');
    imshow(frame);
    title('Adjust sync ROI, then double-click inside the rectangle');
    
    if exist('drawrectangle', 'file') == 2
        roiHandle = drawrectangle('Position', roi, 'Color', 'r');
        wait(roiHandle);
        roi = roiHandle.Position;
    elseif exist('imrect', 'file') == 2
        roiHandle = imrect(gca, roi);
        roi = wait(roiHandle);
    else
        close(figureHandle);
        error('extractSyncFromCamera:RoiToolMissing', ...
            'ROI inspection requires drawrectangle or imrect.');
    end
    
    if ishandle(figureHandle)
        close(figureHandle);
    end
end

function grayFrame = frameToGray(frame)
    if ndims(frame) == 3
        frame = double(frame);
        grayFrame = 0.2989 .* frame(:,:,1) + ...
            0.5870 .* frame(:,:,2) + ...
            0.1140 .* frame(:,:,3);
    else
        grayFrame = double(frame);
    end
    end
    
    function roi = clampRoi(roi, frameSize)
    roi = round(double(roi));
    
    if any(~isfinite(roi)) || roi(3) <= 0 || roi(4) <= 0
        error('extractSyncFromCamera:InvalidRoi', ...
            'ROI must be [x y width height] with positive width and height.');
    end
    
    frameHeight = frameSize(1);
    frameWidth = frameSize(2);
    
    x1 = max(1, roi(1));
    y1 = max(1, roi(2));
    x2 = min(frameWidth, roi(1) + roi(3) - 1);
    y2 = min(frameHeight, roi(2) + roi(4) - 1);
    
    if x1 > x2 || y1 > y2
        error('extractSyncFromCamera:RoiOutOfBounds', ...
            'ROI [%d %d %d %d] does not overlap the video frame.', roi);
    end
    
    roi = [x1, y1, x2 - x1 + 1, y2 - y1 + 1];
end


function quotedPath = shellQuote(pathIn)
    pathIn = char(pathIn);
    pathIn = pathIn(:).';   % force row vector for horzcat

    sq = char(39);          % single quote: '
    dq = char(34);          % double quote: "

    quotedPath = [sq strrep(pathIn, sq, [sq dq sq dq sq]) sq];
end