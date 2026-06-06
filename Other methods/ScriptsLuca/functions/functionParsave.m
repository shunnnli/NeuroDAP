function functionParsaveFigure(fig, plotFullPath)
     
    saveas(fig, plotFullPath, 'pdf')
    saveas(fig, plotFullPath, 'fig')
    
end