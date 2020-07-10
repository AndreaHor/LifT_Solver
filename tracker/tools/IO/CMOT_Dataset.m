classdef (Abstract) CMOT_Dataset < handle
    
    
   properties (Access = protected)
        datasetName_
        datasetSubsetName_
        datasetFolder_
        sequenceProperties_
        configurations_
   end
   
   methods (Abstract)
       folderName = set_dataset_subset(obj,datasetRootDir,datasetSubsetName)
       file_loader = input_file_loader(obj)
   end
   
   
   methods (Access = public)
       
       function mot_sequence = load_sequence(obj,sequenceName)
           sequenceProperties = obj.initSequenceProperties(sequenceName);
           obj.sequenceProperties_ = sequenceProperties;
           mot_sequence = CMOT_Sequence(obj);
           mot_sequence = mot_sequence.load_detections();
       end
       
       function configurations = get_configurations(obj)
           configurations = obj.configurations_;
       end
       
       function sequenceProperties = get_sequence_properties(obj)
           sequenceProperties = obj.sequenceProperties_;
       end

       function sequenceProperties = set_sequence_properties(obj,sequenceProperties)
           obj.sequenceProperties_ = sequenceProperties;
       end
   end
   
   methods (Access = private)
       
       function sequenceProperties = initSequenceProperties(obj,sequenceName)
           sequenceFolder = fullfile(obj.datasetFolder_,sequenceName);
           sequenceProperties = load_ini_as_struct(fullfile(sequenceFolder,'seqinfo.ini'));
           sequenceProperties = sequenceProperties.Sequence;
           sequenceProperties.imDir = fullfile(sequenceFolder,sequenceProperties.imDir,'/');
           sequenceProperties.Dir = sequenceFolder;
       end 
   end
   
   methods (Static, Access = public)
      
       function obj = load_dataset(datasetName,isTrain,configuration)
           
           if contains(datasetName,'MOT1')
              obj = CMOT_MOTChallenge(); 
           else
              error('Dataset not implemented') 
           end
           
           obj.datasetName_= datasetName;
           obj.configurations_ = configuration;
           obj.set_dataset_subset(isTrain);
           
       end
       

       
   end
   
end