function create_folder(folderName)

try
    mkdir(folderName)
    if ~isfolder(folderName)
        mkdir(folderName)
    end
catch
    disp(['error while creating folder ',folderName,' .'])
end
end