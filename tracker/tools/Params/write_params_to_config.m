function write_params_to_config(params,filename)

    ini = IniConfig();
    
    fields = fieldnames(params);
    
    for iField = 1:numel(fields)
       group = fields{iField};
       ini.AddSections(group);
       keys = fieldnames(params.(group));
       for key = 1:numel(keys)
        ini.AddKeys(group,keys{key},params.(group).(keys{key}));   
       end
    end
    ini.WriteFile(filename);
end
