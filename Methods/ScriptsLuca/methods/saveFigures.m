function saveFigures(fig, name, path, options)

arguments
    fig
    name
    path

    options.resolution double = 300
    options.savePNG logical = true
    options.savePDF logical = false
    options.saveFIG logical = false
end

makeFigureFolder = false;

if isempty(dir(fullfile(path))); mkdir(path); end

if sum([options.savePNG, options.savePDF, options.saveFIG]) > 1
    % Make figure folder if more than 1 file type needs to be saved
    mkdir(path,name);
    makeFigureFolder = true;
end

if options.saveFIG
    if ~makeFigureFolder; saveas(fig,fullfile(path,strcat(name,'.fig'))); 
    else; saveas(fig,fullfile(path,name,strcat(name,'.fig'))); 
    end
end
if options.savePNG
    if ~makeFigureFolder; exportgraphics(fig, fullfile(path,strcat(name,'.png')), 'Resolution', options.resolution); 
    else; exportgraphics(fig, fullfile(path,name,strcat(name,'.png')), 'Resolution', options.resolution);
    end
end
if options.savePDF
    if ~makeFigureFolder; exportgraphics(fig, fullfile(path,strcat(name,'.pdf')), 'ContentType', 'vector'); 
    else; exportgraphics(fig, fullfile(path,name,strcat(name,'.pdf')), 'ContentType', 'vector');
    end
end

disp(strcat('Finished: figure "',name,'" saved'));

end