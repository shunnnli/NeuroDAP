function amp = getAmplitude(maxPerTrial,minPerTrial)

arguments
    maxPerTrial double
    minPerTrial double
end

amp = (abs(maxPerTrial)>=abs(minPerTrial)).*maxPerTrial + (abs(maxPerTrial) < abs(minPerTrial)).*minPerTrial;

end