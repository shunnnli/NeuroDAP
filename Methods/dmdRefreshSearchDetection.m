function curCell = dmdRefreshSearchDetection(curCell,searchIdx)
% Re-score cached traces so an old cells_DMD table cannot retain old QC rules.
r = curCell.("Response map"){1};
cm = r.currentMap{searchIdx};
qc = curCell.QC{1}{searchIdx};
if isfield(r,'repMap') && numel(r.repMap)>=searchIdx && ~isempty(r.repMap{searchIdx})
    rm = r.repMap{searchIdx};
else
    rm = cell(size(cm));
    for d = 1:numel(cm)
        counts = cellfun(@(x) size(x,1),cm{d});
        rm{d} = cell(size(cm{d}));
        if isempty(qc{d}) || ~isfield(qc{d},'Sweep') || isempty(qc{d}.Sweep)
            for j=1:numel(counts); rm{d}{j} = "untracked_"+string((1:counts(j))'); end
        elseif all(counts(counts>0)==numel(qc{d}.Sweep))
            % QC rows and trace rows both follow stable acquisition order.
            names = string(qc{d}.Sweep(:));
            for j=1:numel(counts); rm{d}{j} = names(1:counts(j)); end
        else
            % Unequal repetitions require the original per-row sweep identity.
            sessionPath = curCell.Session;
            if ~isfolder(sessionPath); sessionPath = osPathSwitch(sessionPath); end
            path = fullfile(sessionPath,['cell',num2str(curCell.Cell)], ...
                [curCell.Epochs{1}{searchIdx},'_depth',num2str(r.depths{searchIdx}(d)),'.mat']);
            saved = load(path,'spotsAtDepth');
            sequence = r.spotSequence{searchIdx}{d};
            assert(height(saved.spotsAtDepth)==numel(sequence), ...
                'dmdRefreshSearchDetection:rowMismatch','Stored spot sequence does not match source rows.');
            for j=1:numel(counts)
                rm{d}{j} = string(saved.spotsAtDepth.Sweep(sequence==j));
                assert(numel(rm{d}{j})==counts(j), ...
                    'dmdRefreshSearchDetection:countMismatch','Stored sweep identities do not match trace rows.');
            end
        end
    end
end
opt = curCell.Options{1};
criteria = dmdDefaultCriteria();
opt.minSweeps = criteria.minSweeps;
opt.fallbackMaxRs = criteria.fallbackMaxRs;
detection = dmdScoreSearchDetection(cm,rm,qc,opt,opt.timeRange);
for field = {'chargeExcessZ','chargeExcessPc','derivativeExcessZ', ...
        'reanalysisHotspot','nSweepsKept','keptRows','repMap'}
    r.(field{1}){searchIdx} = detection.(field{1});
end
curCell.("Response map"){1} = r;
curCell.Stats{1}.filtered{searchIdx} = detection.filtered;
end
