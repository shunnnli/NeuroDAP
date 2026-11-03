function tests = testSortResultsFolders

tests = functiontests(localfunctions);

end


function testNewestIsLastAcrossAMonthBoundary(testCase)

% dir's own date string is '28-Oct-2026 ...' / '03-Nov-2026 ...', and sorting
% those as text puts the November folder FIRST, so taking {end} as "newest"
% returned the October one. That is how a freshly written Results folder gets
% skipped in favour of a stale one.
root = localMakeFolders(testCase,{'Results-20261028','Results-20261103'}, ...
                                 {'202610281743','202611030912'});

list = sortResultsFolders(root);

verifyEqual(testCase,list{end,1},'Results-20261103');

end


function testOrderIsOldestFirst(testCase)

root = localMakeFolders(testCase,{'Results-20261103','Results-20260723','Results-20261029'}, ...
                                 {'202611030912','202607231200','202610290800'});

list = sortResultsFolders(root);

verifyEqual(testCase,list(:,1), ...
            {'Results-20260723';'Results-20261029';'Results-20261103'});

end


function testMissingParentGivesEmptyList(testCase)

verifyTrue(testCase,isempty(sortResultsFolders(tempname)));

end


function testNonResultsFoldersAreIgnored(testCase)

root = localMakeFolders(testCase,{'Results-20261028','Epochs-20261103'}, ...
                                 {'202610281743','202611030912'});

list = sortResultsFolders(root);

verifyEqual(testCase,size(list,1),1);
verifyEqual(testCase,list{end,1},'Results-20261028');

end


%% Helpers

function root = localMakeFolders(testCase,names,stamps)

% Real folders with real modification times, because the ordering being tested
% is the one dir() reports.
assumeTrue(testCase,isunix,'Needs touch -t to set folder dates.');

root = tempname;
mkdir(root);
addTeardown(testCase,@() rmdir(root,'s'));

for k = 1:numel(names)
    mkdir(fullfile(root,names{k}));
    system(sprintf('touch -t %s "%s"',stamps{k},fullfile(root,names{k})));
end

end
