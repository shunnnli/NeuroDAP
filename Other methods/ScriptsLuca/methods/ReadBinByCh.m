% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [1 X nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.
%
function dataArray = ReadBinByCh(samp0, nSamp, meta, binName, path, ch)

    
    % Channel params
    nChan = str2double(meta.nSavedChans);
    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    samp0 = max(samp0, 0);
    nSamp = floor(min(nSamp, nFileSamp - samp0));
    sizeA = [1, nSamp];
    precision = 'int16=>double';
    skip = 2 * (nChan-1);

    fid = fopen(fullfile(path, binName), 'rb');
    fseek(fid, samp0 * 2 * nChan + 2 * (ch-1), 'bof'); %start reading from the specified channel
    dataArray = fread(fid, sizeA, precision, skip); %skip so that we are not reading other channels
    fclose(fid);


end % ReadBinByCh