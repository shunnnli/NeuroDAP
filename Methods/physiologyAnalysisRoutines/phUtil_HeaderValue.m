function out=phUtil_HeaderValue(header, varName, convert)
	if nargin<3
		convert=0;
	end
	pos=strfind(header, [varName '=']);
	if isempty(pos)
		disp(['WARNING: phUtil_HeaderValue: No header field: ' varName]);
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
		if convert
			out=str2double(out);
			if isempty(out)
				out=nan;
			end
			
		end
	end
		
		
	