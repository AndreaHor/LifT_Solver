function params=ReadParamsFromConfig(ini)

sections = ini.GetSections();
params=struct();

for ii=1:numel(sections)
    [keys, count_keys] = ini.GetKeys(sections{ii});
    values = ini.GetValues(sections{ii}, keys);
    currentSection = sections{ii}(2:end-1);
    params_section = struct();
    for kk=1:count_keys
        params_section=setfield(params_section,keys{kk},values{kk});
    end
    params = setfield(params,currentSection,params_section);
end

end
