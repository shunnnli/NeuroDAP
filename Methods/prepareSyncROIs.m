function prepareSyncROIs(sessionList)
    for s = 1:length(sessionList)
        cameraFiles = dir(fullfile(sessionList{s}, 'cam*.avi'));
        if length(cameraFiles) ~= 1
            error('prepareSyncROIs:CameraFileCount', ...
                'Expected one cam*.avi file in %s, but found %d.', ...
                sessionList{s}, length(cameraFiles));
        end

        fprintf('Select sync ROI for RTPP session: %s\n', sessionList{s});
        cameraPath = fullfile(cameraFiles.folder, cameraFiles.name);
        extractSyncFromCamera(cameraPath, inspectROI=true, selectROIOnly=true);
    end
end
