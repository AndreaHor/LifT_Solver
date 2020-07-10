classdef CMOT_Sequence
    
    properties(Access = public)
        sequenceName_
        detections_
        
        
    end
    
    properties(Access = protected)
        datasetObject_
    end
    
    methods (Access = public)
        
        function obj = CMOT_Sequence(datasetObj)
            obj.datasetObject_ = datasetObj;
            obj.sequenceName_ = obj.get_sequence_properties().name;
        end
        
        function obj = load_detections(obj)
            sequenceProperties = obj.get_sequence_properties();
            csvFile = fullfile(sequenceProperties.Dir,'det','det_preprocessed.csv');
            inputDetectionFile = obj.datasetObject_.input_file_loader(csvFile);
            detections = inputDetectionFile.loadDetections();
            imageEnumerationFormat = inputDetectionFile.get_image_enumeration_format();
            
            detections.descriptor_.identifier = [1:numel(detections.descriptor_.fn)]';
            sequenceProperties.imageEnumerationFormat = imageEnumerationFormat;
            obj.set_sequence_properties(sequenceProperties);
            obj.detections_ = detections;
        end
        
        function configurations = get_configurations(obj)
            configurations = obj.datasetObject_.get_configurations();
        end
        
        function sequenceProperties = get_sequence_properties(obj)
            sequenceProperties = obj.datasetObject_.get_sequence_properties();
        end

        function sequenceProperties = set_sequence_properties(obj,sequenceProperties)
            sequenceProperties = obj.datasetObject_.set_sequence_properties(sequenceProperties);
        end
    end
    
    
    
end