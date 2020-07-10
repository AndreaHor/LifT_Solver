classdef CNodes < handle
    
    
    properties
        descriptor_
        boxFormat
    end
    
    
    methods (Static, Access = private)
        
        function struct = convertLowerMid2UpperLeft(struct)
            x = [struct.ximg];
            y = [struct.yimg];
            w = [struct.width];
            h = [struct.height];
            x = x-0.5*w;
            y = y-h;
            struct.ximg = x;
            struct.yimg = y;
        end
    end

    methods (Access = private)
        function descriptor = getDescriptorFromIdx(obj,idx)
            descriptor = obj.descriptor_;
            for fields = fieldnames(descriptor)'
                currentField = [descriptor.(fields{1})];
                if ~isempty(currentField)
                    if isnumeric(currentField)
                        if numel(currentField) == length(currentField)
                            %currentField is a vector
                            currentField = currentField(idx);
                        elseif strcmp(fields,'emb')
                            currentField = currentField(idx,:);
                        else
                            error('indexing for matrix fields not implemented, for field %s',fields{1})
                        end
                    elseif iscell(currentField)
                        currentField = {currentField{idx}};
                    else
                        error('indexing for field %s not implemented ',fields{1})
                    end
                    descriptor.(fields{1}) = currentField;
                end
            end
        end
        function descriptor = getDescriptorFromFrames(obj,frames)
            [fn,~,~] = obj.frameRange;
            idx = ismember(fn,frames);
            descriptor = obj.getDescriptorFromIdx(idx);
        end
        
        function descriptor = getDescriptorFromBoxID(obj,boxID)
           idx = ismember([obj.descriptor_.identifier],boxID);
           descriptor = obj.getDescriptorFromIdx(idx);
        end
        
        function decriptor = getDescriptorFromLabel(obj,labels)
            
            idx = ismember([obj.descriptor_.label],labels);
            decriptor = obj.getDescriptorFromIdx(idx);
        end
    end
    
    methods (Access = public, Static = true)
       
        
        function saveAsCSV(trajectories,outputFile)
            
            label_all = [];
            x_all = [];
            y_all = [];
            w_all = [];
            h_all = [];
            fn_all = [];
            pdet_all = [];
            for iTraj = 1:numel(trajectories)
                trajectory = trajectories(iTraj);
                [x,y,w,h,fn] = trajectory.getCoordinatesUpperLeft();
                
                label = [trajectory.descriptor_.label];
                pdet = [trajectory.descriptor_.pdet];
                label_all = [label_all; label];
                x_all = [x_all; x];
                y_all = [y_all; y];
                w_all = [w_all; w];
                h_all = [h_all; h];
                fn_all = [fn_all; fn];
                pdet_all = [pdet_all; pdet];
            end
            
            %detMat = [fn,label,x,y,w,h,pdet,x3D,y3D,z3D];
            detMat = [fn_all,label_all,x_all,y_all,w_all,h_all,pdet_all];
            csvwrite([outputFile,'.txt'],detMat);
        end
        
        function tracklet = mergeTrajectories(trajectories)
            
            descriptor = [];
            strFields = fieldnames(trajectories(1).descriptor_);
            
            for field = 1:numel(strFields)
                mergedField = [];
                for iTraj = 1:numel(trajectories)
                    currentField = trajectories(iTraj).descriptor_.(strFields{field});
                    if numel(currentField)>1
                        if size(currentField,1)>size(currentField,2)
                            mergedField = [mergedField; currentField];
                        else
                            mergedField = [mergedField, currentField];
                        end
                    else
                        if size(mergedField,1)>=size(mergedField,2)
                            mergedField =  [mergedField;currentField];
                        else
                            mergedField =  [mergedField, currentField];
                        end
                    end
                end
                descriptor.(strFields{field}) = mergedField;
            end
            [~,sortIdx] = sort(descriptor.fn);
            
            for field = 1:numel(strFields)
                currentField = descriptor.(strFields{field});
               
                descriptor.(strFields{field}) = currentField(sortIdx);
            end
            
            tracklet = CNodes(descriptor,trajectories(1).boxFormat);
            
        end
        
    end
    
    methods (Access = public)
        
        function obj = CNodes(descriptor,boxCoordinatesFormat)
            if strcmp(boxCoordinatesFormat,'upper_left')
                obj.descriptor_ = descriptor;
            elseif strcmp(boxCoordinatesFormat,'lower_mid')
                obj.descriptor_ = CNodes.convertLowerMid2UpperLeft(descriptor);
            end
            obj.boxFormat = 'upper_left';
        end
        function trackletList = splitByLabel(obj)
            labels = obj.descriptor_.label;
            uniqueLabels = unique(labels);
            
            trackletList = [];
            for iLabel = 1:numel(uniqueLabels)
                idxCurrentLabel = find(labels  == uniqueLabels(iLabel));
                if ~isempty(trackletList)
                    trackletList = [trackletList, obj.selectDetectionsFromIdx(idxCurrentLabel)];
                else
                    trackletList = obj.selectDetectionsFromIdx(idxCurrentLabel);
                end
            end
        end

        function [x,y,w,h,fn] = getCoordinatesUpperLeft(obj)
            x = [obj.descriptor_.ximg];
            y = [obj.descriptor_.yimg];
            w = [obj.descriptor_.width];
            h = [obj.descriptor_.height];
            fn = [obj.descriptor_.fn];
        end

        function [labels] = getLabels(obj)
            labels = [obj.descriptor_.label];
        end
        
        function uids = getIdentifiers(obj)
           uids = [obj.descriptor_.identifier];             
        end
        
        function [x,y,w,h,fn] = getCoordinatesLowerMid(obj)
            [x,y,w,h,fn] = obj.getCoordinatesUpperLeft();
            x = x+0.5*w;
            y = y+h;
        end
        
        function [x,y,w,h,fn] =  getCoordinatesCenter(obj)
            [x,y,w,h,fn] = obj.getCoordinatesUpperLeft();
            y = y+0.5*h;
            x = x+0.5*w;
        end
        
        function newTracklet = selectDetectionsFromIdx(obj,idx)
            newTracklet = CNodes(obj.getDescriptorFromIdx(idx),obj.boxFormat);
        end
        
        function newTracklet = selectDetectionsFromBoxID(obj,idx)
            newTracklet = CNodes(obj.getDescriptorFromBoxID(idx),obj.boxFormat);
        end

        
        function [fn,minFn,maxFn] = frameRange(obj)
            fn = [obj.descriptor_.fn];
            minFn = min(fn);
            maxFn = max(fn);
        end
        
        function newDetections = selectDetectionsFromFrames(obj,frames)
            newDetections = CNodes(obj.getDescriptorFromFrames(frames),obj.boxFormat);
        end
        
        function nDetections = length(obj)
            nDetections = numel([obj.descriptor_.ximg]);
        end

        function newDetections = selectDetectionsFromLabel(obj,labels)
            newDetections = CNodes(obj.getDescriptorFromLabel(labels),obj.boxFormat);
        end
 
    end
end
