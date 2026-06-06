function zone=getZone(newCell)
	if nargin<1 || isempty(newCell) % this is a cludge - is empty comes in, return the max number of zones
		%zone=4;
		zone=2;
	else
		% use this for the CTB defined zone
		%zoneLabel=newCell.Injection;
		zoneLabel=newCell.type;
		
		% use this for the anatomically defined zone
		% if newCell.ML<=0
		% 	zoneLabel='M';
		% elseif newCell.ML<=400
		% 	zoneLabel='C';
		% else
		% 	zoneLabel='L';
		% end

		if strcmp(zoneLabel, 'gfp+')
			zone=1;
        elseif strcmp(zoneLabel, 'gfp-')
			zone=2;
	end
	
% 		if strcmp(zoneLabel, 'M')
% 			zone=2;
% 		elseif strcmp(zoneLabel, 'C')
% 			zone=3;
% 		elseif strcmp(zoneLabel, 'L')
% 			zone=4;
% 		elseif strcmp(zoneLabel, 'NA')
% 			zoneLabel=newCell.Notes;
% 
% 			if strcmp(zoneLabel, 'M')
% 				zone=2;
% 			elseif strcmp(zoneLabel, 'C')
% 				zone=3;
% 			elseif strcmp(zoneLabel, 'L')
% 				zone=4;			
% 			end
% 		end
	end