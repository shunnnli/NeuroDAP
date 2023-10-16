function val=phUtil_parsePulsePatternString(pString, param)
    val=[];

    ff=strfind(pString, param);
    if isempty(ff)
        disp([param ' not found in ' pString]);
    else
       ff=ff(1);
       ffEqual=strfind(pString(ff+1:end), '=');
       ffSemi=strfind(pString(ff+1:end), ';');
       if isempty(ffEqual) || isempty(ffSemi)
            disp([param ' malformed in ' pString]);
       else
           vString=pString(ffEqual+ff+1:ffSemi+ff-1);
       end
        val=str2double(vString);
    end
end