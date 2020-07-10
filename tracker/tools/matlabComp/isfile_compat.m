function fileAvail= isfile_compat(dirName)
    fileAvail = exist(dirName) == 2;
end
