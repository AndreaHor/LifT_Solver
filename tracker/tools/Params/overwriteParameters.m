function paramsTracking = overwriteParameters(paramsTracking,overwriteList)    

 for entry = 1:numel(overwriteList)
    paramsTracking.(overwriteList(entry).group).(overwriteList(entry).key) = overwriteList(entry).value;
 end

end
