function params = load_ini_as_struct(filename)
ini_path = IniConfig;
ini_path.ReadFile(filename);
params = ReadParamsFromConfig(ini_path);

end