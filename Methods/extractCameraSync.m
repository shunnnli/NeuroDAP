function extractCameraSync(sessionpath, options)

    arguments
        sessionpath string
        options.sessionName string
        options.fileFormat = '.avi'
    end

    %% Video -> CSV (Sync ROI intensity, FrameCounter, FrameTime)    
    filepath = dir(fullfile(sessionpath,['*',options.fileFormat]));
    if size(filepath,1) > 1
        error('More than 1 camera recording found. Skipped!');
    end

    filename = fullfile(filepath.folder, filepath.name);
    v = VideoReader(filename);
    
    %% --- Select sync ROI on the first frame ---
    firstFrame = readFrame(v);
    
    figure('Name','Select Sync ROI');
    imshow(firstFrame);
    title('Draw rectangle over sync light area, then double-click inside it');
    h = drawrectangle('Color','r');
    wait(h);                        % wait until ROI is finalized
    roi = round(h.Position);         % [x y w h]
    close;
    
    % Rewind to start (important after reading first frame)
    v.CurrentTime = 0;
    
    %% --- Loop through frames and compute outputs ---
    syncVals   = [];   % mean intensity in ROI
    frameCount = [];
    frameTime  = [];
    
    i = 0;
    while hasFrame(v)
        % VideoReader.CurrentTime is the timestamp (seconds) of the NEXT frame to read
        t = v.CurrentTime;
    
        frame = readFrame(v);
        i = i + 1;
    
        % Convert to grayscale if needed
        if size(frame,3) == 3
            gray = rgb2gray(frame);
        else
            gray = frame;
        end
    
        % ROI coordinates (clamped to image bounds)
        x = roi(1); y = roi(2); w = roi(3); hgt = roi(4);
        x1 = max(1, x);
        y1 = max(1, y);
        x2 = min(size(gray,2), x + w - 1);
        y2 = min(size(gray,1), y + hgt - 1);
    
        region = gray(y1:y2, x1:x2);
    
        % Sync intensity: mean pixel value in ROI (use sum(region(:)) if you prefer)
        syncVals(i,1) = mean(double(region(:)));
    
        frameCount(i,1) = i;
        frameTime(i,1)  = t;
    end
    
    %% --- Save to CSV ---
    outCsv   = fullfile(sessionpath, "times_cam1.csv");
    T = table(syncVals, frameCount, frameTime, ...
              'VariableNames', {'Sync','FrameCounter','FrameTime'});
    writetable(T, outCsv);
    
    fprintf("Wrote %d frames to %s\n", height(T), outCsv);

end