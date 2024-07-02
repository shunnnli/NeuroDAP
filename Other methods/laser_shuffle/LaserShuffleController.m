classdef LaserShuffleController
    %LaserShuffleController Main class for dealing with Laser Shuffle UI    
    
    properties
        dataDir; %directory where all the data files are
        dataFiles; %the valid data files in dataDir
        
        fileHeaders; %the header columns in each data file
        fileName2Index; %map from data file name to index in dataFiles
        
        laserName; %name of selected laser
        
        availableCells; %map of cells per file, key is filename, value is array of cell names
        
        parameters; %map of parameters and their values
        
    end
    
    methods
        
        %initializer method
        function obj = LaserShuffleController(obj)

            dparams = containers.Map();
            dparams('laser_name') = 'Laser_1';
            dparams('bin_size') = 5;
            dparams('baseline_start') = -750;
            dparams('baseline_end') = -250;
            dparams('analysis_start') = 0;
            dparams('analysis_end') = 100;
            dparams('pause_between') = 0;
            dparams('sig_latency') = 15;
            
            obj.parameters = dparams;
            
        end
        
        %set the data directory, look for valid files
        function obj = setDataDir(obj, dataDirToUse)
           obj.dataDir = dataDirToUse;
           [dfiles,fheaders] = obj.loadFiles(obj.dataDir);
           obj.dataFiles = dfiles;
           obj.fileHeaders = fheaders;
           obj.fileName2Index = containers.Map();
           for k = 1:length(dfiles)              
               obj.fileName2Index(dfiles{k}) = k;               
           end
           
           if isempty(obj.dataFiles) 
               errordlg('The data directory has no data files!');              
           end         
        end
        
        %check to see if column name is valid
        function filesWithCol = filesWithColumn(obj, colName)
            
            filesWithCol = {};
            for k = 1:length(obj.fileHeaders)
                for j = 1:length(obj.fileHeaders{k})
                    %fprintf('comparing %s with %s\n', obj.fileHeaders{k}{j}, colName)
                    if strcmp(obj.fileHeaders{k}{j}, colName)                       
                        filesWithCol{end+1} = obj.dataFiles{k};                        
                    end
                end
            end
        end
               
        %set the Laser name and find available cells
        function obj = setLaserName(obj, lname)
        
            filesWithCol = obj.filesWithColumn(lname);
    
            if isempty(filesWithCol)
                errordlg(sprintf('Laser name %s not found in any of the files!', lname));
                obj.laserName = '';
                return;        
            end
            
            obj.laserName = lname;
            filesToSearch = obj.filesWithColumn(lname);
            
            cnames = containers.Map();
            for k = 1:length(filesToSearch)
                
                fname = filesToSearch{k};
                findex = obj.fileName2Index(fname);
                fheaders = obj.fileHeaders{findex};
                
                cells = LaserShuffleController.getCellNames(fheaders, lname);
                if ~isempty(cells)
                    cnames(fname) = cells;
                end
            end            
            
            obj.availableCells = cnames;
            
        end
        
        %run the laser shuffle analysis
        function obj = runAnalysis(obj, selectedFiles, selectedCellsPerFile)
            
            paramVals.numRandShifts = 500; % Number of times to shuffle spike times from each trial
            paramVals.percentileCutOff = 0.995; % Cut-off for determining a specific bin is significant
            paramVals.randShiftMultiplier = 1; % (sec) Determines the magnitude of time randomized time shift
            paramVals.pauseDur = 0; % (sec) Duration to pause after each plot; 0=no pause; -1 = wait for user
            paramVals.minSpikesPerSigBin = 3;% 0.005; % Discard cells that don't have enough spikes to judge analysis
            paramVals.minSigBins = 3; % Minimum number of bins that must be significant in order to count cell
            paramVals.histHistBins = 0:1:5000; % Set bin size and range for histogram of PSTH values
            paramVals.numDropFractions = 50; % Number of times to recompute the score omitting a fraction of the pulses
            paramVals.selectCells = 0; % 1=query user for filename & cell; 0=analyze all files & cells
            paramVals.highlightBinTimes = [-0.1 0.005 0.015 0.1]; % Bin times to highligh on plots
                        
            paramVals.pauseDur = obj.parameters('pause_between');
            paramVals.selectCells = {};
            paramVals.laserColName = obj.laserName;
            paramVals.psthWindow = [-1 2];
            paramVals.psthBinSize = obj.parameters('bin_size')*1e-3;
            paramVals.preLaserWindow = [obj.parameters('baseline_start') obj.parameters('baseline_end')]*1e-3;
            paramVals.laserPulseWindow = [obj.parameters('analysis_start') obj.parameters('analysis_end')]*1e-3;
            
            paramVals.psthBins = paramVals.psthWindow(1):paramVals.psthBinSize:paramVals.psthWindow(2);
            paramVals.laserPSTHBins = find(paramVals.psthBins >= paramVals.laserPulseWindow(1) & paramVals.psthBins < paramVals.laserPulseWindow(2));
            paramVals.laserPSTHBins = (paramVals.laserPSTHBins(1)-1):paramVals.laserPSTHBins(end);

            allOutVals = {};
            for k = 1:length(selectedFiles)               
                selectedCells = selectedCellsPerFile{k};
                fullFileName = fullfile(obj.dataDir, selectedFiles{k});                
                outVals = laserShuffleAnalysis(fullFileName, paramVals, selectedCells);
                allOutVals{end+1} = outVals;
            end
            
            %write the output values to a file
            outDir = fullfile(obj.dataDir, 'output');
            if ~exist(outDir, 'dir')
                %try to make output directory
                [success, msg, msgid] = mkdir(outDir);
                if ~success
                    errdlg(sprintf('Could not create output directory at %s, defaulting to %s', ...
                           outDir, obj.dataDir));
                    outDir = obj.dataDir;
                end
            end
            
            d = datestr(clock());
            d = strrep(d, ' ', '_');
            d = strrep(d, ':', '.');
            fileName = sprintf('output_%s.csv', d);
            outFile = fullfile(outDir, fileName);            
            f = fopen(outFile, 'w');
            
            %write parameter values
            fprintf(f, 'Laser,BinSize,BaselineStart,BaselineEnd,AnalysisStart,AnalysisEnd,SigLatency\n');
            fprintf(f, '%s,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n', obj.laserName, obj.parameters('bin_size'), ...
                    obj.parameters('baseline_start'), obj.parameters('baseline_end'), ...
                    obj.parameters('analysis_start'), obj.parameters('analysis_end'), ...
                    obj.parameters('sig_latency'));
            fprintf(f, '\n');
            
            %write output value headers
            fprintf(f, 'FileName,CellName,FirstSigBin,LaserModulated\n');
            
            %write output values
            for k = 1:length(allOutVals)
                outVals = allOutVals{k};                
                for j = 1:outVals.numCells
                    fsmb = outVals.firstSigMaxBin(j);
                    isSignificant = (fsmb > 0) & (fsmb <= (obj.parameters('sig_latency')*1e-3 + 1e-3));
                    
                    fprintf(f, '%s,%s,%0.3f,%d\n', ...
                            selectedFiles{k}, outVals.cellNames{j}, ...
                            fsmb, isSignificant);                    
                end                
            end            
            fclose(f);
            display(sprintf('Wrote output file to %s', outFile));
            
        end
        
    end
    
    methods (Static)
        
        %get the cells from a file with a specified laser name
        function cNames = getCellNames(fHeaders, lName)            
            cNames = {};
            passed = 0;
            for k = 1:length(fHeaders)
                if passed
                    cNames{end+1} = fHeaders{k};
                end                
                if strcmp(fHeaders{k}, lName)
                    passed = 1;
                end
            end            
        end
        
        function [goodFiles,fheaders] = loadFiles(dataDir)
            fprintf('Loading files from %s\n', dataDir);
            fnames = dir(dataDir);
            goodFiles = {};
            fheaders = {};
            for k = 1:length(fnames)
                if ~fnames(k).isdir
                    fileName = fnames(k).name;
                    fullFileName = fullfile(dataDir, fileName);
                    headers = LaserShuffleController.getFileHeaders(fullFileName);
                    if ~isempty(headers) > 0
                        goodFiles{end+1} = fileName;
                        fheaders{end+1} = headers;
                    end                    
                end                
            end                     
        end
        
        %read the first line of a file to get the headers
        function headerNames = getFileHeaders(fileName)
            headerNames = {};
            f = fopen(fileName, 'r');
            firstLine = fgetl(f);
            if ~isempty(firstLine)
                cols = textscan(firstLine, '%s', 'delimiter', ',');
                if length(cols{1}) > 1
                    headerNames = cols{1};
                end                
            end            
            fclose(f);
        end

    end
    
end

