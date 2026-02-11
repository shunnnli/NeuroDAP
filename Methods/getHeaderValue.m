function out = getHeaderValue(header, varName, options)

% Shun Li, 7/1/2024
% Modified from phUtil_HeaderString by Bernardo Sabatini

arguments
    header
    varName
    
    options.pulseVar string

    options.convert logical = false
    options.removeQuote logical = true
end


if ~isempty(varName)
    % Extract value
    pos=strfind(header, [varName '=']);

    if isempty(pos)
	    disp(['WARNING: getHeaderValue: No header field: ' varName]);
	    out=[];
	    if convert
		    out=nane;
	    end
    else
	    posEq=strfind(header(pos:end), '=');
	    posRt=strfind(header(pos:end), 13);
        if isempty(posRt)
	        out=header(pos+posEq(1):end);
	    else
		    out=header(pos+posEq(1):pos+posRt(1)-2);
        end
    
        % Remove quote if necessary (i.e. in cycleName)
        if options.removeQuote && out(1) == "'" && out(end) == "'"
            out = out(2:end-1);
        end
    
        % Convert output to double
	    if options.convert
		    out=str2double(out);
		    if isempty(out); out=nan; end
	    end
    end
end

% Extract pulse value if necessary
if isfield(options,'pulseVar')
    if isempty(varName); out = header; end
    pString = out; out = [];
    ff = strfind(pString, options.pulseVar);

    if isempty(ff)
        disp(strcat(options.pulseVar, ' not found in ', pString));
    else
       ff = ff(1);
       ffEqual = strfind(pString(ff+1:end), '=');
       ffSemi = strfind(pString(ff+1:end), ';');
       if isempty(ffEqual) || isempty(ffSemi)
            disp(strcat(param,' malformed in ',pString));
       else
           startIndex = ffEqual + ff + 1; startIndex = startIndex(1);
           endIndex = ffSemi + ff - 1; endIndex = endIndex(1);
           vString = pString(startIndex:endIndex);
       end
        out = str2double(vString);
    end
end

end