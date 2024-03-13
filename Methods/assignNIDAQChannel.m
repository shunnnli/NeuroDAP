function channels = assignNIDAQChannel(analogNI, digitalNI, options)

arguments
    analogNI double
    digitalNI double
    options.user string = 'Shun'
    options.inverStim logical = true
    options.getConsecutive logical = true
end

if strcmpi('Shun',options.user)
    % Convert other analog signals to digital 
    % Rising edge: temp == 1
    % Falling edge: temp == -1
    leftLick = (analogNI(1,:)> (max(analogNI(1,:))/2));
    temp = [false, diff(leftLick)];
    leftLick = (temp==-1);
    
    rightLick = (analogNI(2,:)> (max(analogNI(2,:))/2));
    temp = [false, diff(rightLick)];
    rightLick = (temp==-1);
    
    blink = (analogNI(3,:)> (max(analogNI(3,:))/2));
    temp = [false, diff(blink)];
    blinkRaw = (temp==1);
    blink = find(blinkRaw);
    
    gyro = single(analogNI(4:6,:));
    photometry_raw = analogNI(8,:); % Photometry PMT
    
    % Digital channels
    if options.invertStim
        digitalNI(7,:) = ~digitalNI(7,:);
        digitalNI(8,:) = ~digitalNI(8,:);
    end

    temp = false(size(digitalNI));
    for i=1:size(digitalNI,1)
        temp2 = [false, diff(digitalNI(i,:))];
        temp(i,:) = (temp2==1);
    end
    
    % Calculate length of each pulse
    if options.getConsecutive
        airpuff = temp(1,:) .* getConsecutive(digitalNI(1,:))./nidq.Fs;
        leftSolenoid = temp(3,:) .* getConsecutive(digitalNI(3,:))./nidq.Fs;
        rightSolenoid = temp(4,:) .* getConsecutive(digitalNI(4,:))./nidq.Fs;
        leftTone = temp(5,:) .* getConsecutive(digitalNI(5,:))./nidq.Fs;
        rightTone = temp(6,:) .* getConsecutive(digitalNI(6,:))./nidq.Fs;
        redLaser = temp(7,:) .* getConsecutive(digitalNI(7,:))./nidq.Fs;
        blueLaser = temp(8,:) .* getConsecutive(digitalNI(8,:))./nidq.Fs;
    else
        airpuff = temp(1,:);
        leftSolenoid = temp(3,:);
        rightSolenoid = temp(4,:);
        leftTone = temp(5,:);
        rightTone = temp(6,:);
        redLaser = temp(7,:);
        blueLaser = temp(8,:);
    end

    % Set up allTrials: not including omissionTrials
    % 1 is left trial, 2 is right trial
    allTones = zeros(size(leftTone));
    allTones(leftTone ~= 0) = 1;
    allTones(rightTone ~= 0) = 2;

    syncNI = digitalNI(2,:); % Save sync channel from NIDAQ data separately

    % Set output
    channels{1} = leftLick; channels{2} = rightLick;
    channels{3} = blink; channels{4} = gyro; channels{5} = photometry_raw;

    channels{6} = airpuff; channels{7} = leftSolenoid; channels{8} = rightSolenoid;
    channels{9} = leftTone; channels{10} = rightTone; channels{11} = redLaser; 
    channels{12} = blueLaser; channels{13} = allTones;
    channels{14} = syncNI;
    return

elseif strcmpi('Kevin',options.user)
    % Convert other analog signals to digital 
    % Rising edge: temp == 1
    % Falling edge: temp == -1

    syncNI = (analogNI(2,:) > (max(analogNI(2,:))/2));
    temp = [false, diff(syncNI)];
    syncNI = (temp==-1);

    leftLick = (analogNI(3,:)> (max(analogNI(3,:))/2));
    temp = [false, diff(leftLick)];
    leftLick = (temp==-1);
    
    rightLick = (analogNI(4,:)> (max(analogNI(4,:))/2));
    temp = [false, diff(rightLick)];
    rightLick = (temp==-1);
    
    LeftNose = (analogNI(5,:)> (max(analogNI(5,:))/2));
    temp = [false, diff(LeftNose)];
    LeftNose = (temp==1);
    
    % Digital channels
    temp = false(size(digitalNI));
    for i=1:size(digitalNI,1)
        temp2 = [false, diff(digitalNI(i,:))];
        temp(i,:) = (temp2==1);
    end
    
    % Calculate length of each pulse
    rightLED = temp(1,:);
    rightSolenoid = temp(2,:);
    rightLick = temp(3,:);
    rightNose = temp(4,:);
    centerLED = temp(5,:);
    centerNose = temp(6,:);
    leftLED = temp(7,:);
    leftSolenoid = temp(8,:);

    % Set output
    channels{1} = syncNI; channels{2} = leftLick;
    channels{3} = rightLick; channels{4} = LeftNose;

    channels{5} = rightLED; channels{6} = rightSolenoid; channels{7} = rightLick;
    channels{8} = rightNose; channels{9} = centerLED; channels{10} = centerNose; 
    channels{11} = leftLED; channels{12} = leftSolenoid;
    return

end