classdef LaserShuffleView
    %LaserShuffleView takes care of all the UI stuff for Laser Shuffle
    %   Detailed explanation goes here
    
    properties        
        cellFiles; %cell array for the file associated with each cell in the listbox
    end
    
    methods
        
        function obj = setDataFiles(obj, controller, handles)            
            set(handles.listbox_valid_files, 'String', controller.dataFiles);
            set(handles.text_data_dir_value, 'String', controller.dataDir);
        end
        
        function sfiles = getSelectedFiles(obj, controller, handles)           
            %get the selected files
            sfiles_index = get(handles.listbox_valid_files, 'Value');
            validfiles = get(handles.listbox_valid_files, 'String');
            sfiles = cell(length(sfiles_index), 1);
            for k = 1:length(sfiles_index)
                sfiles{k} = validfiles{sfiles_index(k)};
            end            
        end
        
        function obj = selectAllCells(obj, controller, handles)
            
            svals = get(handles.listbox_cells, 'String');
            if ~isempty(svals)
                svals_index = 1:length(svals);
                set(handles.listbox_cells, 'Value', svals_index);
            end
        end
        
        function obj = selectAllFiles(obj, controller, handles)
            
            svals = get(handles.listbox_valid_files, 'String');
            if ~isempty(svals)
                svals_index = 1:length(svals);
                set(handles.listbox_valid_files, 'Value', svals_index);
            end
        end
        
        function obj = setAvailableCells(obj, controller, handles)
            
            if isempty(controller.availableCells)
                return;
            end
            
            sfiles = obj.getSelectedFiles(controller, handles);
            
            %check the available files against the selected files and build
            %up a list of cells
            allCells = {};            
            cfiles = {};
            for k = 1:length(sfiles)
                fkey = sfiles{k};                
                if controller.availableCells.isKey(fkey)
                    clist = controller.availableCells(fkey);
                    for j = 1:length(clist)                       
                        allCells{end+1} = clist{j};
                        cfiles{end+1} = fkey;
                    end
                end                
            end
            obj.cellFiles = cfiles;
            
            %set the listbox with the values
            set(handles.listbox_cells, 'String', allCells);            
        end
        
        function cellsPerFile = getSelectedCellsPerFile(obj, controller, handles)
            
            sfiles = obj.getSelectedFiles(controller, handles);
            
            fileMap = containers.Map(); %key is filename, value is list of cells
            selIndices = get(handles.listbox_cells, 'Value');
            cnames = get(handles.listbox_cells, 'String');
            
            for k = 1:length(selIndices)
                selIndex = selIndices(k);
                cName = cnames{selIndex};
                cFile = obj.cellFiles{selIndex};
                if ~fileMap.isKey(cFile)
                    fileMap(cFile) = {};
                end
                ca = fileMap(cFile);
                ca{end+1} = cName;
                fileMap(cFile) = ca;
            end
            
            cellsPerFile = cell(length(sfiles), 1);
            for k = 1:length(cellsPerFile)
                fname = sfiles{k};
                cellsPerFile{k} = fileMap(fname);
            end            
        end
        
        function obj = setDefaults(obj, controller, handles)
            
            p = controller.parameters;
            
            set(handles.edit_laser_name, 'String', p('laser_name'));
            
            set(handles.edit_bin_size, 'Value', p('bin_size'));
            set(handles.edit_bin_size, 'String', sprintf('%0.1f', p('bin_size')));
            
            set(handles.edit_baseline_start, 'Value', p('baseline_start'));
            set(handles.edit_baseline_start, 'String', sprintf('%0.1f', p('baseline_start')));
            
            set(handles.edit_baseline_end, 'Value', p('baseline_end'));
            set(handles.edit_baseline_end, 'String', sprintf('%0.1f', p('baseline_end')));
            
            set(handles.edit_analysis_start, 'Value', p('analysis_start'));
            set(handles.edit_analysis_start, 'String', sprintf('%0.1f', p('analysis_start')));
            
            set(handles.edit_analysis_end, 'Value', p('analysis_end'));
            set(handles.edit_analysis_end, 'String', sprintf('%0.1f', p('analysis_end')));
            
            set(handles.edit_sig_latency, 'Value', p('sig_latency'));
            set(handles.edit_sig_latency, 'String', sprintf('%0.1f', p('sig_latency')));
                                    
        end        
    end
    
end

