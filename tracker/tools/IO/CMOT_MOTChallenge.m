classdef CMOT_MOTChallenge < CMOT_Dataset
   
    methods (Access = public)
        
        function datasetSubsetName = set_dataset_subset(obj,isTrain)
            if isTrain
                datasetSubsetName = 'train';
            else
                datasetSubsetName = 'test';
            end
            obj.datasetFolder_ = fullfile(obj.configurations_.Data.datasetsDir,obj.datasetName_,datasetSubsetName,'/');
            
        end
        
        function file_loader = input_file_loader(obj,filename)
            file_loader = CMOTFile(filename);
        end

    end
end