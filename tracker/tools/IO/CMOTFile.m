classdef CMOTFile < CTrackingFile
    % Tracking file in MOTChallenge format
    
    properties
       type2column = {'fn','label','ximg','yimg','width','height','pdet','x3d','y3d','z3d'};
    end
    
    methods(Access = public)
        
        function obj = CMOTFile(filename)
            obj.setFileName(filename);
        end
        
        function tracklets = loadDetections(obj)
            tracklets = obj.loadFromTable();
        end
        
        
        function boxCoordinatesFormat = get_box_coordinates_format(obj)
            boxCoordinatesFormat = 'upper_left';
        end
        
        function imEnumFormat = get_image_enumeration_format(obj)
            imEnumFormat = '%06d';
        end
    

        function idx = columnIdx(obj,columnSelector)
            idx = find(strcmp(obj.type2column,columnSelector));
            if isempty(idx)
                error('Column %s is not part of MOT file format.',columnSelector);
            end
        end
    end
end

