classdef CVisualizer < handle
    
    properties  (Access = public, Constant = true)
        colorizer = distinguishable_colors(30);
    end
    
    properties (Access = public, Constant = false)
        imEnumFormat_
        imExt_
        imPath_
        imCounter_
        saveIm_
        saveDir_
        figHandle_
        
        traceInfo_
        trace_
    end
    
    methods (Access=public, Static = true)
        
        function img = plotText(img,text,position)
            if nargin == 2
                position = [10 40];
            end
            img = insertText(img,position,text);
        end
        
        function img = plotFrameNumber(img,fn,position)
            if nargin == 2
                position = [10 10];
            end
            img = insertText(img,position,['Frame = ',num2str(fn)]);
        end
        
        function img = plotFrameNumberFromDetection(img,detection,position)
            assert(detection.length == 1,'Function can be called only for single detections');
            [fn,~,~] = detection.frameRange;
            if nargin == 2
                img = CVisualizer.plotFrameNumber(img,fn);
            else
                img = CVisualizer.plotFrameNumber(img,fn,position);
            end
        end
        
        function img = plotDetections(img,detections,isLabeled,text)
            
            [x,y,w,h,~] = detections.getCoordinatesUpperLeft;
            
            if nargin <4
                if isLabeled
                    ids = detections.getLabels();
                    idModulo = mod(ids,30)+1;
                    color = round(255*CVisualizer.colorizer(idModulo,:));
                    img = insertObjectAnnotation(img,'Rectangle',[x,y,w,h],ids,'LineWidth',5,'Color',color);
                else
                    img  = insertShape(img,'Rectangle',[x,y,w,h],'LineWidth',5);
                end
            else
                if isLabeled
                    ids = detections.getLabels();
                    idModulo = mod(ids,30)+1;
                    color = round(255*CVisualizer.colorizer(idModulo,:));
                    img = insertObjectAnnotation(img,'Rectangle',[x,y,w,h],text,'LineWidth',5,'Color',color);
                else
                    img  = insertShape(img,'Rectangle',[x,y,w,h],'LineWidth',5);
                end
            end
        end
        
        function img = plotDetectionsAsCircle(img,detections,isLabeled,text)
            
            [x,y,w,h,~] = detections.getCoordinatesUpperLeft();
            img  = insertShape(img,'FilledRectangle',[x,y,w,h],'LineWidth',5,'Opacity',0.3);
        end
        
    end
    
    methods (Access = public)
        
        
        function newImg = loadFrame(obj,fr)
            newImg = imread(fullfile(obj.imPath_,[num2str(fr,obj.imEnumFormat_),obj.imExt_]));
        end
        
        function newImg = loadFrameOfDetection(obj,detection)
            [fr,~,~] = detection.frameRange;
            assert(numel(fr) == 1,'single detection expected')
            newImg= obj.loadFrame(fr);
        end
        
        function drawImg(obj,img)
            figHandle = obj.figHandle_;
            activeHandle= ~isempty(figHandle) && ishandle(figHandle) && strcmp(get(figHandle,'type'),'image');
            if ~activeHandle & ~obj.saveIm_
                figure;
                obj.figHandle_ = imshow(img);
                hold on
            else
                obj.figHandle_.CData = img;
            end
            
            if obj.saveIm_
                obj.imCounter_ = obj.imCounter_+1;
                imwrite(img,[obj.saveDir_,num2str(obj.imCounter_,obj.imEnumFormat_),'.png']);
            else
                pause(0.001)
            end
        end
        
        function img = drawTrace(obj,img,detections)
            [~,fn,~] = detections.frameRange;
            [x,y,~,~,~] = detections.getCoordinatesLowerMid;
            ids = detections.getLabels();
            idModulo = mod(ids,30)+1;
            color = round(255*CVisualizer.colorizer(idModulo,:));
            trace = obj.trace_;
            traceInfo = obj.traceInfo_;
            
            positions = [];
            colorAll = [];
            activeTrajs = [];
            for iID = 1:numel(ids)
                activeTrajs = [activeTrajs, ids(iID)];
                idxCurrentDet = ids == ids(iID);
                if numel(trace)>= ids(iID)
                    currentTrace = trace{ids(iID)};
                    currentTraceInfo = traceInfo{ids(iID)};
                else
                    currentTrace = [];
                    currentTraceInfo = [];
                end
                rmRange = currentTraceInfo <= fn-30;
                currentTrace(rmRange,:) = [];
                currentTraceInfo(rmRange) = [];
                currentTrace = [currentTrace; x(idxCurrentDet),y(idxCurrentDet)];
                positions = [positions; currentTrace];
                colorAll = [colorAll; repmat(color(iID,:),size(currentTrace,1),1)];
                currentTraceInfo = [currentTraceInfo, fn];
                trace{ids(iID)} = currentTrace;
                traceInfo{ids(iID)} = currentTraceInfo;
            end
            obj.trace_ = trace;
            obj.traceInfo_ = traceInfo;
            positions(:,3) = 4;
            
            for iTr = 1:numel(activeTrajs)
                pos_lin = trace{activeTrajs(iTr)};
                pos_lin = pos_lin';
                pos_lin = pos_lin(:);
                pos_lin = pos_lin';
                if numel(pos_lin)>2
                    img = insertShape(img,'Line',pos_lin,'Color',color(iTr,:),'LineWidth',3);
                else
                    img = insertShape(img,'Circle',[pos_lin, 4],'Color',color(iTr,:));
                end
            end
            %
        end
    end
    
    methods (Access = public)
        function obj = CVisualizer(mot_sequence,saveImg,saveDir)
            sequenceProperties = mot_sequence.get_sequence_properties();
            if saveImg
                obj.saveIm_ = true;
                obj.saveDir_ = saveDir;
            else
                obj.saveIm_ = false;
            end
            obj.imCounter_= 0;
            obj.imExt_ = sequenceProperties.imExt;
            obj.imEnumFormat_ = sequenceProperties.imageEnumerationFormat;
            obj.imPath_ = sequenceProperties.imDir;
            obj.trace_ = [];
            obj.traceInfo_ = [];
        end
        
        function reset(obj)
            obj.imCounter_ = 0;
            obj.figHandle_ = [];
        end
        
        function drawTrajectories(obj,trajectories,detections)
            mergedTrajectories = CNodes.mergeTrajectories(trajectories);
            [~,minFN,maxFn] = mergedTrajectories.frameRange;
            for fr = minFN:maxFn
                detectionsInFrame = mergedTrajectories.selectDetectionsFromFrames(fr);
                if ~isempty(detections)
                    inputdetectionsInFrame = detections.selectDetectionsFromFrames(fr);
                end
                
                
                img = obj.loadFrame(fr);
                img = CVisualizer.plotFrameNumber(img,fr);
                if detectionsInFrame.length>0
                    img = CVisualizer.plotDetections(img,detectionsInFrame,true);
                    img = obj.drawTrace(img,detectionsInFrame);
                end
                
                if ~isempty(detections)
                    for idx =1:inputdetectionsInFrame.length
                        input_det_current = inputdetectionsInFrame.selectDetectionsFromIdx(idx);
                        [x_i,y_i,w_i,h_i,~] = input_det_current.getCoordinatesLowerMid;
                        max_iou= 0;
                        for j_idx = 1:detectionsInFrame.length
                            det_tracker = detectionsInFrame.selectDetectionsFromIdx(j_idx);
                            [x_d,y_d,w_d,h_d,~] = det_tracker.getCoordinatesLowerMid;
                            
                            iou_bb = bboxOverlapRatio([x_i,y_i,w_i,h_i],[x_d,y_d,w_d,h_d]);
                            max_iou = max([iou_bb,max_iou]);
                        end
                        if max_iou < 0.8
                            img = CVisualizer.plotDetectionsAsCircle(img,input_det_current,false);
                        end
                    end
                end
                obj.drawImg(img);
                pause(0.05);
            end
        end
        
        function drawTracklet(obj,tracklet)
            nDetections = tracklet.length;
            for iDet = 1:nDetections
                obj.drawDetection(tracklet.selectDetectionsFromIdx(iDet),false);
            end
        end
        
        
        function drawDetectionOnFrame(obj,detection,isLabeled,text,frame)
        
            assert(detection.length == 1,'This method can be called only for a single detection');
            img = obj.loadFrame(frame);
            img = CVisualizer.plotFrameNumberFromDetection(img,detection);
            if nargin >=4
                    img = CVisualizer.plotText(img,text);
            end
            img = CVisualizer.plotDetections(img,detection,isLabeled);
            obj.drawImg(img);
            
        end
        function drawDetection(obj,detection,isLabeled,text,position)
            assert(detection.length == 1,'This method can be called only for a single detection');
            img = obj.loadFrameOfDetection(detection);
            img = CVisualizer.plotFrameNumberFromDetection(img,detection);
            if nargin >=4
                if nargin == 5
                    img = CVisualizer.plotText(img,text,position);
                else
                    img = CVisualizer.plotText(img,text);
                end
            end
            img = CVisualizer.plotDetections(img,detection,isLabeled);
            obj.drawImg(img);
        end
        
    end
    
end
