function sdata = readLaserShuffleData(filePath, laserColumnName)

    [allData, delim] = importdata(filePath);
    ncols = length(allData.colheaders);
    
    %% find laser column index
    laserCol = -1;
    for k = 1:ncols
        if strcmp(allData.colheaders{k}, laserColumnName)
           laserCol = k; 
        end               
    end
    
    if laserCol == -1
        error('Cannot find column for laser name %s in file %s', laserColumnName, filePath);
    end
    
    %% get laser times
    laserTimes = allData.data(:, laserCol);
    laserTimes = laserTimes(~isnan(laserTimes));
    
    %% get cell names and indices
    cell2column = containers.Map();        
    for k = (laserCol+1):ncols        
        cname = allData.colheaders{k};
        cell2column(cname) = k;                
    end    
    
    %% get cell data    
    cellNames = cell2column.keys();
    cellSpikes = cell(length(cellNames), 1);
    for k = 1:length(cellNames)       
        cname = cellNames{k};
        ccol = cell2column(cname);        
        cdata = allData.data(:, ccol);
        cdata = cdata(~isnan(cdata));
        cellSpikes{k} = cdata;
    end
    
    %% clean up
    %clear allData;
    
    %% build output structure
    sdata = struct;
    sdata.laser = laserColumnName;
    sdata.laserTimes = laserTimes;    
    sdata.cellNames = cellNames;
    sdata.cellSpikes = cellSpikes;
    sdata.allData = allData;
    sdata.cell2column = cell2column;
        
end
