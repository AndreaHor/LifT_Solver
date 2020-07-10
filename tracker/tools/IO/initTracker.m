function [params,setup] = initTracker(pathToConfigFiles)

    paramsName = 'params.ini';
    configName = 'configuration.ini';
    
    params = load_ini_as_struct(fullfile(pathToConfigFiles,paramsName));
    setup = load_ini_as_struct(fullfile(pathToConfigFiles,configName));
end
