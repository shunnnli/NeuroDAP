function stats = getCellStats(combined_cells,current_direction,options)

arguments
    combined_cells table
    current_direction string % "exci" or "inhi"
    options.type string = 'auc'

    options.animalRange string = 'All'
    options.cellList double % row in combined_cells
    
    options.included logical = true

    options.average logical = false
    options.withQC logical = false
end

% Define cell to plot
if ~isfield(options,'cellList')
    cellList = 1:size(combined_cells,1);
    if ~strcmpi(options.animalRange,'All')
        cellList = cellList(contains(combined_cells.Animal,options.animalRange));
    end
else
    cellList = options.cellList;
end


for c = 1:length(cellList)
    cellIdx = cellList(c);
    curCell = combined_cells(cellIdx,:); 

    % Get rows
    if contains(current_direction,'exci')
        vholdRows = curCell.Vhold{1} == -70;
    elseif contains(current_direction,'inhi')
        vholdRows = curCell.Vhold{1} >= 0;
    end

    % Get included
    included = vertcat(curCell.Included{1}{vholdRows});

    % Get stats
    if ~strcmpi(options.type,'qc')
        baseline_stats = vertcat(curCell.Stats{1}.baseline.(options.type){vholdRows});
        response_stats = vertcat(curCell.Stats{1}.response.(options.type){vholdRows});
        
        % Filter by included and average if necessary
        if options.included
            baseline_stats = baseline_stats(included == 1);
            response_stats = response_stats(included == 1);
        end
        if options.average
            baseline_stats = mean(baseline_stats);
            response_stats = mean(response_stats);
        end
    else
        qc = mergeStructs(curCell.QC{1}(vholdRows));
        metrics = fieldnames(qc);
        if options.included
            for f = 1:numel(metrics)
                if strcmpi(metrics{f},'included')
                    continue
                end

                data = qc.(metrics{f});

                if length(data) == length(included)
                    data = data(included == 1);
                end
                if options.average
                    data = mean(data,'omitmissing');
                end
                
                qc.(metrics{f}) = data;
            end
        end
    end

    % Store in stats struct
    stats(c).Animal = curCell.Animal;
    stats(c).Task = curCell.Task;
    stats(c).Cell = curCell.Cell;
    stats(c).Vhold = curCell.Vhold{1}(vholdRows);
    if strcmpi(options.type,'qc')
        stats(c).QC = qc;
    else
        stats(c).baseline = baseline_stats;
        stats(c).response = response_stats;
    end
end

end