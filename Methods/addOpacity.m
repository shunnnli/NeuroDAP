function newColor = addOpacity(color,opacity)

if iscell(color)
    newColor = color;
    for i = 1:length(color)
        newColor{i} = 1 - opacity*(1-color{i});
    end
else
    newColor = 1 - opacity*(1-color);
end

end