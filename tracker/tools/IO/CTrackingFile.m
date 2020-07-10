classdef (Abstract) CTrackingFile < handle
    
    properties
        filename_;
    end
        
    methods (Abstract)
        idx = columnIdx(obj,columnSelector)
        
        boxCoordinatesFormat   = get_box_coordinates_format(obj);
        
        imEnumFormat = get_image_enumeration_format(obj);
        
    end

    methods (Access = protected)
        
        
        function detection_nodes =loadFromTable(obj)
            %for input detections with IoU Score and sorted
            filename = obj.getFileName();
            
            detections = readtable(filename,'ReadVariableNames',1,'Delimiter',',');
            
            boxCoordinatesFormat = obj.get_box_coordinates_format();
            
            %inputs are unlabeled detections
            descriptor.fn   = detections.Frame;
            descriptor.ximg = detections.X_UL;
            descriptor.yimg = detections.Y_UL;
            descriptor.width = detections.W;
            descriptor.height = detections.H;
            descriptor.pdet = detections.PDET;
            descriptor.label = -1*ones(size(detections.PDET));
            
            detection_nodes = CNodes(descriptor,boxCoordinatesFormat);
        end
        
        function setFileName(obj,fileName)
            if ~isfile_compat(fileName)
                error('File %s could not be found',fileName);
            end
            obj.filename_ = fileName;
        end
        
        function filename = getFileName(obj)
            if ~isempty(obj.filename_)
                filename = obj.filename_;
            else
                error('Filename not set')
            end
        end
    end
end
