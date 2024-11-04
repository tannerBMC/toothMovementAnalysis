% program for registering surface data for 
% the aligner performance case study for two patients reported in 
% 
% I. Filippon, C. Tanner, J. A. von Jackowski, G. Schulz, T. Toepper, B.
% Mueller, University of Basel, 2024
% Aligner-induced tooth movements in three dimensions using clinical data
% of two patients, Oral, 2024
%
% Weekly oral scans for first nine weeks and scan at end of treatment
% n=1:11, where n=1  is Week 0 before treatment
%               n=10 is Week 9 and 
%               n=11 is at the end of treatment
% 
% Semiautomatic segmentation of teeth in end state (endN=10=refN)
% Manual definition of occlusion plane in end state (endN=10=refN)
%
% Performs rigid and then non-rigid ICP registration
% rigid:      global alignment of all teeth
% non-rigid:  local alignment of individual teeth
%
% evaluate displacements from state n->endN->simN on teeth
%                             
% based on  non-rigid motion, i.e. non-rigid versus rigid ICP results 
% defined only for moving mesh hence need to deform reference 
% and displacements components are then in moving space
%
% run0: do n=[1:refN-1 refN+1:N]   endN=refN
%       propagate semiautomatic teeth segmentation and occlusion plane
%       from state endN=refN to all states n=[1:refN-1 refN+1:N] 
%       via surface registration
%       
% run1: do n=1 endN=[2:N]
%       determine motion from state n=1 to endN=2:N
%       and difference between endN and simN
%       via surface registration
%
% freecad was used to look and cut stl files
%
% May 2024, Christine Tanner, Biomaterials Science Center, University of Basel

clear all
close all

% define paths for other programs
path(path,['utils' filesep])

saveFigs = 0;            % 1: save figures
saveMainFigs = 0;        % 1: save only most important figures for manuscript
                         %    Fig2*png is Fig. 4 in manuscript etc.
saveGAbstractFigs = 0;   % 1: save figures for graphical abstract
makeClean = 0;           % 1: create clean figures for manuscript
 
saveData = 1;            % 1: save teeth errors, achieved, planned, missing
       
buccLowerToRight = 1;    % 1: in Fig5 ensure buccular points to right for lower jaw
doXFlipLowerDisp = 1;    % 1: flip X in output png files for lower jaw
                                    
redoMovie = 0;           % 1: redo movie
       
runThrough = 1;          % 1: run through without breakpoints
            
%% load data definitions and registration parameters
dataDefinitionsCaseStudy

% define/create registration directories
regFigsDir = [figsDir regParStr filesep];
if ~exist(regFigsDir,'file')
    mkdir(regFigsDir);
end

resultDir = [baseDir 'results' filesep regParStr filesep];
if ~exist(resultDir,'file')
    mkdir(resultDir);
end

%% loop over state n->endN and simState n->endN
%
% (run0) do n=[1:refN-1 refN+1:N]   endN=refN, e.g. n=9 endN=refN
% (run1) do n=1 to endN=2:N                    e.g. n=1 endN=9
%
% needed??
% (run2) do n=[2:refN-1 refN+1] to end=refN:refN

for p=1:1  % SELECT patient 1: 3485 (patient A) 2: 6457 (patient B)
    
for n=1:1  %[11 2:9];  %1:N-1 %2:N-1; %N-1;  % SELECT state [1:9 11], i.e. not end state N=refN

endN = 9; %3:2:10 %2:N;  % SELECT end state, first do endN=refN to propagate segmentations

doEndSimRegV = 0:1;  % 0: state_n 1: simState_endN

% do always End state as fixed and register stateN and simStateN to it
if endN==refN
    endStr='End';
else
    endStr=['End' num2str(endN)];
end

% compare same state of planning
simEndN = endN;

%%
for doUpper = 0:0   % 0: lower, 1: upper jaw
    
    for doEndSimReg = doEndSimRegV  
        %close all
        
        % set seed to make registration reproducible
        rng(10*p+doUpper)
        
        extraStr = '';
        if doSmallerRigidTol
            teethRigidStr = '_TolS';
        else
            teethRigidStr = '';
        end
        if doUpper
            partStr = 'upper';
            % collect filenames
            % selected state n
            stlFileStateName = [baseDir pat{p}.dir pat{p}.stateUpper{n} extraStr '.stl'];
            % propagated teeth segmentation
            if doEndSimReg
                % select simulated state at endN
                stateTeethSegDir = [baseDir pat{p}.simDir pat{p}.stateUpperTeeth{endN}(1:end-7) 'teethSeg_' nrRegParStr teethRigidStr filesep];
            else
                % selected measured state n
                stateTeethSegDir = [baseDir pat{p}.dir pat{p}.stateUpperTeeth{n}(1:end-7) 'teethSeg_' nrRegParStr teethRigidStr filesep];
            end
            stlFileStateNameTeeth = [stateTeethSegDir  pat{p}.stateUpperTeeth{n}(end-6:end) extraStr '.stl'];
            stlFileStateNameToothStr = [stateTeethSegDir pat{p}.stateUpperToothStr{n}(end-5:end) extraStr];
            if ~exist(stateTeethSegDir,'file')
                mkdir(stateTeethSegDir)
            end
            occPlaneStateName = [stlFileStateNameTeeth(1:end-4) '_occPlane.mat'];
            % measured endN state
            stlFileEndName = [baseDir pat{p}.dir pat{p}.stateUpper{endN} extraStr '.stl'];
            if endN==refN
                stlFileEndNameTeeth = [baseDir pat{p}.dir pat{p}.endUpperTeeth extraStr '.stl'];
                stlFileEndNameToothStr = [baseDir pat{p}.dir pat{p}.endUpperToothStr extraStr];
            else
                % need transferring!
                endTeethSegDir = [baseDir pat{p}.dir pat{p}.stateUpperTeeth{endN}(1:end-7) 'teethSeg_' nrRegParStr teethRigidStr  filesep];
                stlFileEndNameTeeth = [endTeethSegDir  pat{p}.stateUpperTeeth{endN}(end-6:end) extraStr '.stl'];
                stlFileEndNameToothStr = [endTeethSegDir pat{p}.stateUpperToothStr{endN}(end-5:end) extraStr];
                if exist(stlFileEndNameTeeth,'file')
                    disp(['use transferred teeth segmentation ' stlFileEndNameTeeth])
                else
                    disp([stlFileEndNameTeeth ' missing!'])
                    disp(['transfer teeth segmentation using n=' num2str(endN) ' and endN=' num2str(refN)])
                    return
                end
            end
            occPlaneEndName = [stlFileEndNameTeeth(1:end-4) '_occPlane.mat'];
            % selected simulated state simEndN
            stlFileSimStateName = [baseDir pat{p}.simDir pat{p}.stateSimUpper{simEndN} '.stl'];
            
            % rotation for upper, change from angle to rad
            rotXMrad=pat{p}.rotXM(:,:,1)/180*pi;  
            rotYMrad=pat{p}.rotYM(:,:,1)/180*pi;  
            rotZMrad=pat{p}.rotZM(:,:,1)/180*pi;  
            
            % define which contour to cut
            if doEndSimReg
                stateTRcut = stateSimUpperCutM(simEndN,p);
                stateTRkeepHigh = stateSimUpperKeepHighM(simEndN,p);
                stateTRStraight = stateSimUpperStraightM(simEndN,p);
                if stateTRcut==-1
                    stlOrigFileSimStateName = [baseDir pat{p}.simDir pat{p}.stateSimUpperOrig{simEndN} extraStr '.stl'];
                end
            else
                stateTRcut = stateUpperCutM(n,p);
                stateTRkeepHigh = stateUpperKeepHighM(n,p);
                stateTRStraight = stateUpperStraightM(n,p);
                if stateTRcut==-1
                    stlOrigFileStateName = [baseDir pat{p}.dir pat{p}.stateUpperOrig{n} extraStr '.stl'];
                end
            end
            endTRcut = stateUpperCutM(endN,p);
            endTRkeepHigh = stateUpperKeepHighM(endN,p);
            if endTRcut==-1
                stlOrigFileEndName = [baseDir pat{p}.dir pat{p}.stateUpperOrig{endN} extraStr '.stl'];
            end
            endTRStraight = stateUpperStraightM(endN,p);
                    
            if endN==refN
                occPlane = upperOccPlane{p};
            else
                % need transferring via non-rigid registration!
                if exist(occPlaneEndName,'file')
                    disp(['use transferred occlusion plane ' occPlaneEndName])
                else
                    disp([occPlaneEndName ' missing!'])
                    disp(['transfer occlusion plane using n=' num2str(endN) ' and endN=' num2str(refN)])
                    return
                end
                load(occPlaneEndName);  % load 'occPlaneState'
                occPlane=occPlaneState;
                if ~runThrough
                    keyboard
                end
            end
        else
            partStr = 'lower';
            % selected measured state n
            stlFileStateName = [baseDir pat{p}.dir pat{p}.stateLower{n} extraStr '.stl'];
            % propagated teeth segmentation, here it is simState
            if doEndSimReg
                % simulation at endN
                stateTeethSegDir = [baseDir pat{p}.simDir pat{p}.stateLowerTeeth{endN}(1:end-7) 'teethSeg_' nrRegParStr teethRigidStr filesep];
            else
                % selected measured state n
                stateTeethSegDir = [baseDir pat{p}.dir pat{p}.stateLowerTeeth{n}(1:end-7) 'teethSeg_' nrRegParStr teethRigidStr filesep];
            end
            stlFileStateNameTeeth = [stateTeethSegDir  pat{p}.stateLowerTeeth{n}(end-6:end) extraStr '.stl'];
            stlFileStateNameToothStr = [stateTeethSegDir pat{p}.stateLowerToothStr{n}(end-5:end) extraStr];
            if ~exist(stateTeethSegDir,'file')
                mkdir(stateTeethSegDir)
            end
            occPlaneStateName = [stlFileStateNameTeeth(1:end-4) '_occPlane.mat'];
            
            % measured end state
            stlFileEndName = [baseDir pat{p}.dir pat{p}.stateLower{endN} extraStr '.stl'];
            if endN==refN
                stlFileEndNameTeeth = [baseDir pat{p}.dir pat{p}.endLowerTeeth extraStr '.stl'];
                stlFileEndNameToothStr = [baseDir pat{p}.dir pat{p}.endLowerToothStr extraStr ];
            else
                endTeethSegDir = [baseDir pat{p}.dir pat{p}.stateLowerTeeth{endN}(1:end-7) 'teethSeg_' nrRegParStr teethRigidStr filesep];
                stlFileEndNameTeeth = [endTeethSegDir  pat{p}.stateLowerTeeth{endN}(end-6:end) extraStr '.stl'];
                stlFileEndNameToothStr = [endTeethSegDir pat{p}.stateLowerToothStr{endN}(end-5:end) extraStr];
                if exist(stlFileEndNameTeeth,'file')
                    disp(['use transferred teeth segmentations, e.g. ' stlFileEndNameTeeth])
                else
                    disp([stlFileEndNameTeeth ' missing!'])
                    disp(['transfer teeth segmentations using n=' num2str(endN) ' and endN=' num2str(refN)])
                    return
                end
            end
            occPlaneEndName = [stlFileEndNameTeeth(1:end-4) '_occPlane.mat'];
            % selected simulated state simEndN
            stlFileSimStateName = [baseDir pat{p}.simDir pat{p}.stateSimLower{simEndN} '.stl'];
            
            % rotation for lower, change from angle to rad
            rotXMrad=pat{p}.rotXM(:,:,2)/180*pi;  
            rotYMrad=pat{p}.rotYM(:,:,2)/180*pi;  
            rotZMrad=pat{p}.rotZM(:,:,2)/180*pi;  
            
            % define which contour to cut
            if doEndSimReg
                stateTRcut = stateSimLowerCutM(simEndN,p);
                stateTRkeepHigh = stateSimLowerKeepHighM(simEndN,p);
                stateTRStraight = stateSimLowerStraightM(simEndN,p);
                if stateTRcut==-1
                    stlOrigFileSimStateName = [baseDir pat{p}.simDir pat{p}.stateSimLowerOrig{simEndN} extraStr '.stl'];
                end
            else
                stateTRcut = stateLowerCutM(n,p);
                stateTRkeepHigh = stateLowerKeepHighM(n,p);
                stateTRStraight = stateLowerStraightM(n,p);
                if stateTRcut==-1
                    stlOrigFileStateName = [baseDir pat{p}.dir pat{p}.stateLowerOrig{n} extraStr '.stl'];
                end
            end
            endTRcut = stateLowerCutM(endN,p);
            endTRkeepHigh = stateLowerKeepHighM(refN,p);
            if endTRcut==-1
                stlOrigFileEndName = [baseDir pat{p}.dir pat{p}.stateLowerOrig{endN} extraStr '.stl'];
            end
            endTRStraight = stateLowerStraightM(endN,p);
            if endN==refN
                occPlane = lowerOccPlane{p};
            else
                if exist(occPlaneEndName,'file')
                    disp(['use transferred occlusion plane ' occPlaneEndName])
                else
                    disp([occPlaneEndName ' missing!'])
                    disp(['transfer occlusion plane using n=' num2str(endN) ' and endN=' num2str(refN)])
                    return
                end
                load(occPlaneEndName);  % load 'occPlaneState'
                occPlane=occPlaneState;
            end
        end
        
        %% read state n stl file, i.e. other state be registered to endN
        if doEndSimReg
            % this is a simulated state
            % find difference between planned and actual end state
            stateStr = ['Sim' num2str(simEndN)];
            
            % TR triangulation with Points and ConnectivityList
            stateTR = stlread(stlFileSimStateName);
            disp(['read stateTR ' stlFileSimStateName])
            
            % ensure same orientation
            theta=rotXMrad(simEndN,2);   % simulated state n
            thetaY=rotYMrad(simEndN,2);  
            thetaZ=rotZMrad(simEndN,2);  
            disp(['simulated state ' num2str(simEndN) ' X ' num2str(theta) ' Y ' num2str(thetaY) ' Z ' num2str(thetaZ)])
            if theta~=0  || thetaZ~=0 || thetaY~=0
                Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
                Ry = [cos(thetaY) 0 sin(thetaY); 0 1 0; -sin(thetaY) 0 cos(thetaY)];
                Rz = [cos(thetaZ) -sin(thetaZ) 0; sin(thetaZ) cos(thetaZ) 0; 0 0 1];
                rotM=Rz*Ry*Rx;
                rotPoints=stateTR.Points*rotM;
                rotTR=triangulation(stateTR.ConnectivityList,rotPoints);
                stateTR=rotTR;
            end
            % extract normals
            pointNormalsStart = vertexNormal(stateTR);
            triNormalsStart = faceNormal(stateTR);
        else
            % measured state n
            stateStr = ['State' num2str(n)];
            
            % TR triangulation with Points and ConnectivityList
            stateTR = stlread(stlFileStateName);
            disp(['read stateTR ' stlFileStateName]);
            if stateTRcut==-1
                % read original mesh
                stateTRorig = stlread(stlOrigFileStateName);
            end
            
            doRotTesting = 0;
            theta=rotXMrad(n,1);  % measured state n
            thetaY=rotYMrad(n,1);
            thetaZ=rotZMrad(n,1);
            disp(['measured state' num2str(n) ' X ' num2str(theta) ' Y ' num2str(thetaY) ' Z ' num2str(thetaZ)])
            if theta~=0  || thetaZ~=0 || thetaY~=0
                Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
                Ry = [cos(thetaY) 0 sin(thetaY); 0 1 0; -sin(thetaY) 0 cos(thetaY)];
                Rz = [cos(thetaZ) -sin(thetaZ) 0; sin(thetaZ) cos(thetaZ) 0; 0 0 1];
                rotM=Rz*Ry*Rx; 
                if doRotTesting
                    % show individual rotation to ease setting them
                    figure(302)
                    for ri=1:4
                        if ri==1
                            rotPoints=stateTR.Points;
                            titleStr = 'orig';
                            viewV=[1 0 0];
                        elseif ri==2
                            rotPoints=stateTR.Points*Rz;
                            titleStr = 'Rz';
                            viewV=[0 0 1];
                        elseif ri==3
                            rotPoints=stateTR.Points*Rz*Ry;
                            titleStr = 'Rz*Ry';
                            viewV=[0 1 0];
                        elseif ri==4
                            rotPoints=stateTR.Points*Rz*Ry*Rx;
                            titleStr = 'Rz*Ry*Rx';
                            viewV=[1 0 0];
                        end
                        rotTR=triangulation(stateTR.ConnectivityList,rotPoints);
                        tmpTR=rotTR;
                        subplot(2,2,ri)
                        fh=trisurf(tmpTR);
                        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                        axis equal
                        axis tight
                        xlabel('x')
                        ylabel('y')
                        zlabel('z')
                        title(titleStr)
                        lightangle(gca,90,30)
                        view(viewV)
                    end
                end
                rotPoints=stateTR.Points*rotM;
                rotTR=triangulation(stateTR.ConnectivityList,rotPoints);
                stateTR=rotTR;
            end
        end
        disp(['***** Processing ' partStr ' ' stateStr ' to ' endStr])
        if doUpper==0
            % fix y-direction to be the same as for upper
            flipYPoints=stateTR.Points;
            flipYPoints(:,2)=-flipYPoints(:,2);
            stateTR=triangulation(stateTR.ConnectivityList,flipYPoints);
            if stateTRcut==-1 && ~doEndSimReg
                % original mesh
                flipYPoints=stateTRorig.Points;
                flipYPoints(:,2)=-flipYPoints(:,2);
                stateTRorig=triangulation(stateTRorig.ConnectivityList,flipYPoints);
            end
        end
        
        %% read end stl file
        
        % TR triangulation with Points and ConnectivityList
        endTR = stlread(stlFileEndName);
        disp(['read endTR ' stlFileEndName]);
        if endTRcut==-1
            endTRorig = stlread(stlOrigFileEndName);
        end
        if doUpper
            for to=1:length(seg{p}.upperToothStrV)
                toothFname = [stlFileEndNameToothStr seg{p}.upperToothStrV{to} '.stl'];
                seg{p}.endToothTR{to} = stlread(toothFname);
                seg{p}.toothStrV{to} = seg{p}.upperToothStrV{to};
                % assign name and color
                tidx=str2num(seg{p}.upperToothStrV{to}(end));
                seg{p}.name{to}=toothCat{tidx}.name;
                seg{p}.short{to}=toothCat{tidx}.short;
                seg{p}.color{to}=toothCat{tidx}.color;
            end
            if endN==refN
                for to=1:length(seg{p}.extraUpperToothStrV)
                    toothFname = [stlFileEndNameToothStr seg{p}.extraUpperToothStrV{to} '.stl'];
                    seg{p}.extraEndToothTR{to} = stlread(toothFname);
                end
                seg{p}.extraToothStrV = seg{p}.extraUpperToothStrV;
            end
        else
            for to=1:length(seg{p}.lowerToothStrV)
                toothFname = [stlFileEndNameToothStr seg{p}.lowerToothStrV{to} '.stl'];
                seg{p}.endToothTR{to} = stlread(toothFname);
                seg{p}.toothStrV{to} = seg{p}.lowerToothStrV{to};
                tidx=str2num(seg{p}.lowerToothStrV{to}(end));
                seg{p}.name{to}=toothCat{tidx}.name;
                seg{p}.short{to}=toothCat{tidx}.short;
                seg{p}.color{to}=toothCat{tidx}.color;
            end
            if endN==refN
                for to=1:length(seg{p}.extraLowerToothStrV)
                    toothFname = [stlFileEndNameToothStr seg{p}.extraLowerToothStrV{to} '.stl'];
                    seg{p}.extraEndToothTR{to} = stlread(toothFname);
                end
                seg{p}.extraToothStrV = seg{p}.extraLowerToothStrV;
            end
        end
        try
            endTeethTR = stlread(stlFileEndNameTeeth);
            disp(['read endTeethTR ' stlFileEndNameTeeth])
        catch
            % not existing yet, compose from segmented teeth
            allPoints = [];
            allConn = [];
            prevNoP = 0;
            for i=1:length(seg{p}.endToothTR)
                allPoints = [allPoints; seg{p}.endToothTR{i}.Points];
                % need to add point number offset
                allConn = [allConn; seg{p}.endToothTR{i}.ConnectivityList+prevNoP];
                prevNoP=prevNoP+size(seg{p}.endToothTR{i}.Points,1);
            end
            % save combined
            endTeethTR = triangulation(allConn,allPoints);
            figure(11)
            fh=trisurf(endTeethTR);
            axis equal
            keyboard
            stlwrite(endTeethTR,stlFileEndNameTeeth,'text');
        end
        if endN==refN
            % also read teeth segmentation propagated to state
            if exist(stlFileStateNameTeeth,'file') 
                stateTeethTR = stlread(stlFileStateNameTeeth);
                disp(['read stateTeethTR ' stlFileStateNameTeeth])
            end
        end
        
        % ensure same orientation for endTeeth
        theta=rotXMrad(endN,1); % measured end state
        thetaY=rotYMrad(endN,1);
        thetaZ=rotZMrad(endN,1);
        disp(['measured end state X ' num2str(theta) ' Y ' num2str(thetaY) ' Z ' num2str(thetaZ)])
        if theta~=0  || thetaZ~=0 || thetaY~=0
            Rx = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            Ry = [cos(thetaY) 0 sin(thetaY); 0 1 0; -sin(thetaY) 0 cos(thetaY)];
            Rz = [cos(thetaZ) -sin(thetaZ) 0; sin(thetaZ) cos(thetaZ) 0; 0 0 1];
            rotM=Rz*Ry*Rx; 
            rotPoints=endTR.Points*rotM;
            rotTR=triangulation(endTR.ConnectivityList,rotPoints);
            endTR=rotTR;
            if endN==refN
                rotPoints=endTeethTR.Points*rotM;
                rotTR=triangulation(endTeethTR.ConnectivityList,rotPoints);
                endTeethTR=rotTR;
                for to=1:length(seg{p}.toothStrV)
                    rotPoints=seg{p}.endToothTR{to}.Points*rotM;
                    rotTR=triangulation(seg{p}.endToothTR{to}.ConnectivityList,rotPoints);
                    seg{p}.endToothTR{to}=rotTR;
                end
                if endN==refN
                    for to=1:length(seg{p}.extraToothStrV)
                        rotPoints=seg{p}.extraEndToothTR{to}.Points*rotM;
                        rotTR=triangulation(seg{p}.extraEndToothTR{to}.ConnectivityList,rotPoints);
                        seg{p}.extraEndToothTR{to}=rotTR;
                    end
                end
            end
            if endTRcut==-1
                rotPoints=endTRorig.Points*rotM;
                rotTR=triangulation(endTRorig.ConnectivityList,rotPoints);
                endTRorig=rotTR;
            end
        end
        
        if doUpper==0 
            % fix y-direction to be the same as for upper
            flipYPoints=endTR.Points;
            flipYPoints(:,2)=-flipYPoints(:,2);
            endTR=triangulation(endTR.ConnectivityList,flipYPoints);
            if endTRcut==-1
                flipYPoints=endTRorig.Points;
                flipYPoints(:,2)=-flipYPoints(:,2);
                endTRorig=triangulation(endTRorig.ConnectivityList,flipYPoints);
            end
            % no flipping of propagated teeth segmentation needed 
            if endN==refN
                flipYPoints=endTeethTR.Points;
                flipYPoints(:,2)=-flipYPoints(:,2);
                endTeethTR=triangulation(endTeethTR.ConnectivityList,flipYPoints);
                for to=1:length(seg{p}.toothStrV)
                    flipYPoints=endTeethTR.Points;
                    flipYPoints(:,2)=-flipYPoints(:,2);
                    endTeethTR=triangulation(endTeethTR.ConnectivityList,flipYPoints);
                    
                    flipYPoints=seg{p}.endToothTR{to}.Points;
                    flipYPoints(:,2)=-flipYPoints(:,2);
                    seg{p}.endToothTR{to}=triangulation(seg{p}.endToothTR{to}.ConnectivityList,flipYPoints);
                end
                if endN==refN
                    for to=1:length(seg{p}.extraToothStrV)
                        flipYPoints=seg{p}.extraEndToothTR{to}.Points;
                        flipYPoints(:,2)=-flipYPoints(:,2);
                        seg{p}.extraEndToothTR{to}=triangulation(seg{p}.extraEndToothTR{to}.ConnectivityList,flipYPoints);
                    end
                end
            end
        end
        
        %% show teeth segmentation for checking it
        
        clear legendStr
        disp('Figure 330: Show teeth segmentation')
        figure(330)
        fh=trisurf(endTR);
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        legendStr{1}='Teeth';
        hold on
        showIdxV=[1:length(seg{p}.endToothTR)]; 
        for i1=1:length(showIdxV)
            i=showIdxV(i1);
            fh=trisurf(seg{p}.endToothTR{i});
            set(fh,'FaceColor',seg{p}.color{i},'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            legendStr{i1+1}=[seg{p}.toothStrV{i} ' ' seg{p}.short{i}];
        end
        hold off
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        if makeClean
            [caz,cel]=view([0 -1 0]);
        else
            [caz,cel]=view([1 -1 0]);
            title([pat{p}.id ' ' endStr ' ' partStr])
            legend(legendStr,'Location','bestoutside')
        end
        lightangle(gca,caz,cel)
        material(materialV)
        if saveFigs
            if makeClean
                axis off
                print(330,'-dpng',[figsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_clean.png'])
            else
                print(330,'-dpng',[figsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_withWisdom.png'])
            end
        end
        
        if ~runThrough
            keyboard
        end
        
        %% define/show 3 points for occlusion plane
        disp('Figure  63: Define/show occlusion plane')
        
        figure(63)
        fh=trisurf(endTR);
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        axis tight
        if doUpper
            xlabel('x (right P2 (1-7) ->left P3 (2-7))')
        else
            xlabel('x (left P3 (3-7) ->right P2 (4-7))')
        end
        ylabel('y')
        zlabel('z')
        if makeClean
            [caz,cel]=view([0 -1 0]);
        else
            [caz,cel]=view([1 -1 0]);
            title([pat{p}.id ' ' endStr ' ' partStr])
        end
        lightangle(gca,caz,cel)
        if isempty(occPlane.P1)
            ax = gca;
            ax.Interactions = [dataTipInteraction rotateInteraction];
            disp('Rotate by holding right mouse down')
            disp('Pick points by right mouse click, box with voxel coordinates X, Y, Z, appears')
            disp('Pick P1 (between teeth), P2 (right side, upper 1-7, lower 4-7), P3 (left side, upper 2-7, lower 3-7)')
            disp('Right click mouse, use Export Cursor Data to Workspace, save in dataDefinitionIgnacio')
            disp('cursor_info.Position')
            keyboard
        else
            hold on
            plot3(occPlane.P1(1),occPlane.P1(2),occPlane.P1(3),'r*','MarkerSize',ms/2,'LineWidth',lw/2)
            plot3(occPlane.P2(1),occPlane.P2(2),occPlane.P2(3),'g*','MarkerSize',ms/2,'LineWidth',lw/2)
            plot3(occPlane.P3(1),occPlane.P3(2),occPlane.P3(3),'b*','MarkerSize',ms/2,'LineWidth',lw/2)
            occPlane.midP = (occPlane.P2+occPlane.P3)/2;
            plot3(occPlane.midP(1),occPlane.midP(2),occPlane.midP(3),'k*','MarkerSize',ms/2,'LineWidth',lw/2)
            line([occPlane.P2(1) occPlane.P3(1)],[occPlane.P2(2) occPlane.P3(2)],[occPlane.P2(3) occPlane.P3(3)],'LineWidth',lw/2)
            line([occPlane.P1(1) occPlane.midP(1)],[occPlane.P1(2) occPlane.midP(2)],[occPlane.P1(3) occPlane.midP(3)],'LineWidth',lw/2)
            hold off
            if doUpper
                legend({'surface','P1','P2 (1-7)','P3 (2-7)','midP'})
            else
                legend({'surface','P1','P2 (4-7)','P3 (3-7)','midP'})
            end
            set(gca,'FontSize',fs)
            if saveFigs
                if makeClean
                    axis off
                    legend off
                    figFname = [figsDir pat{p}.id '_' partStr '_' endStr '_occlusionPlane_clean.png'];
                else
                    figFname = [figsDir pat{p}.id '_' partStr '_' endStr '_occlusionPlane.png'];
                end
                print(63,'-dpng',figFname)
            end
        end
        if ~runThrough
            keyboard
        end
        %% show teeth before cutting
        if ~doEndSimReg
            disp('Figure 301: Show meshes before cutting')
            figure(301)
            subplot(1,2,2)
            fh=trisurf(endTeethTR);
            title([endStr 'Teeth input'])
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            axis equal
            axis tight
            xlabel('x')
            ylabel('y')
            zlabel('z')
            lightangle(gca,90,30)
            view([1 0 0])
            subplot(1,2,1)
            if stateTRcut==-1
                % original mesh
                fh=trisurf(stateTRorig);
            else
                fh=trisurf(stateTR);
            end
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            axis equal
            axis tight
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title([stateStr ' input'])
            lightangle(gca,90,30)
            view([1 0 0])
            if saveFigs
                print(301,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr '_input.png'])
            end
        end
        %% cut mesh, either bei contour line or manually
        
        if doEndSimReg
            keepAmount=10;
        else
            % do not cut as not closed
            keepAmount=0;
        end
        
        if keepAmount>0
            % based on residual amount
            if endTRcut>-1
                endTR = cutMesh(endTR, doUpper, keepAmount, endTRkeepHigh);
            end
            if stateTRcut>-1
                stateTR = cutMesh(stateTR, doUpper, keepAmount, stateTRkeepHigh);
            end
        else
            % based on contour index
            if endTRcut>-1
                endTR = cutMeshIndex(endTR,endTRcut,endTRkeepHigh,endTRStraight);
            end
            if stateTRcut>-1
                stateTR = cutMeshIndex(stateTR,stateTRcut,stateTRkeepHigh,stateTRStraight);
            end
        end
        
        %% show meshes after cut
        disp('Figure   3: Show meshes after cutting')
        figure(3)
        if doEndSimReg
            subplot(1,2,1)
        else
            subplot(1,2,2)
        end
        if doEndSimReg
            fh=trisurf(endTeethTR);
            title([endStr 'Teeth cut'])
        end
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        axis tight
        xlabel('x')
        ylabel('y')
        zlabel('z')
        lightangle(gca,90,30)
        view([1 0 0])
        
        figure(3)
        if doEndSimReg
            subplot(1,2,2)
        else
            subplot(1,2,1)
        end
        fh=trisurf(stateTR);
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        axis tight
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title([stateStr ' cut'])
        lightangle(gca,90,30)
        view([1 0 0])
        if saveFigs
            print(3,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr '_cut.png'])
        end
        
        if ~runThrough
            keyboard
        end
        
        %% make endTR orient with occlusion plane
        if alignOccPlane
            disp('Figure 722: Make endTR orient with occlusion plane')
            if isempty(occPlane.P1)
                disp('occlusion plane NOT defined!')
                keyboard
            else
                % collect points
                X=[occPlane.P1; occPlane.P2; occPlane.P3; occPlane.midP];
                X3P=[occPlane.P1; occPlane.P2; occPlane.P3];
                Xo=mean(X3P);              % mean position of 3 points
                
                % determine normal to occlusion plane
                % this should become the y vector
                ptCloudPlane = pointCloud(X3P);
                modelPlane = pcfitplane(ptCloudPlane,0.1);
                if modelPlane.Normal(2)<0
                    yV=-modelPlane.Normal;
                else
                    yV=modelPlane.Normal;
                end
                
                % line P2 to P3 should become x vector
                xV=(occPlane.P2-occPlane.P3);
                xV=xV/norm(xV);
                if xV(1)<0
                    xV=-xV;
                end
                % determine orthogonal x-vector
                zV=cross(yV,xV);
                if zV(3)<0
                    zV=-zV;
                end
                
                % rotation matrix to align with these axis
                newR = [xV' yV' zV'];
                
                % apply rotation to points
                dX2=bsxfun(@minus,X,Xo);   % position same cente
                newX=dX2*newR+Xo;
                
                % update occlusion plane points
                occPlane.P1 = newX(1,:);
                occPlane.P2 = newX(2,:);
                occPlane.P3 = newX(3,:);
                occPlane.midP = newX(4,:);
                
                % apply rotation to endTR
                dX2=bsxfun(@minus,endTR.Points,Xo);   % position same cente
                rotPoints=dX2*newR+Xo;
                rotTR=triangulation(endTR.ConnectivityList,rotPoints);
                
                figure(722)
                fh=trisurf(endTR);
                set(fh,'FaceColor',BMCblue,'EdgeColor','none','FaceAlpha',fAlpha*0.5,'FaceLighting',fLighting)
                hold on
                fh=trisurf(rotTR);
                set(fh,'FaceColor',BMCred,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                % show occlusion plane lines
                line([occPlane.P2(1) occPlane.P3(1)],[occPlane.P2(2) occPlane.P3(2)],[occPlane.P2(3) occPlane.P3(3)],'color','k')
                line([occPlane.P1(1) occPlane.midP(1)],[occPlane.P1(2) occPlane.midP(2)],[occPlane.P1(3) occPlane.midP(3)],'color','k')
                hold off
                axis equal
                % stick to fixed image axis
                xlabel('x')
                ylabel('y')
                zlabel('z')
                [caz,cel]=view([1 -1 0]);
                lightangle(gca,caz,cel)
                material(materialV)
                legend('orig','alignOccPlane')
                
                % update end state
                endTR=rotTR;
                % update segmented teeth
                if ~isempty(seg{p}.toothStrV)
                    X2=endTeethTR.Points;       % mesh points
                    dX2=bsxfun(@minus,X2,Xo);   % position same cente
                    Y2=dX2*newR+Xo;
                    endTeethTR = triangulation(endTeethTR.ConnectivityList,Y2);
                    for to=1:length(seg{p}.toothStrV)
                        X2=seg{p}.endToothTR{to}.Points;  % mesh points
                        dX2=bsxfun(@minus,X2,Xo);         % position same cente
                        Y2=dX2*newR+Xo;
                        seg{p}.endToothTR{to} = triangulation(seg{p}.endToothTR{to}.ConnectivityList,Y2);
                    end
                    % 
                    if endN==refN
                        for to=1:length(seg{p}.extraToothStrV)
                            X2=seg{p}.extraEndToothTR{to}.Points;  % mesh points
                            dX2=bsxfun(@minus,X2,Xo);         % position same cente
                            Y2=dX2*newR+Xo;
                            seg{p}.extraEndToothTR{to} = triangulation(seg{p}.extraEndToothTR{to}.ConnectivityList,Y2);
                        end
                    end
                    
                end
            end
        end
        if ~runThrough
            keyboard
        end
        %% registration of surfaces
        disp('Registration of surfaces')
        
        % rigid ICP data
        fixCloudFull = pointCloud(endTR.Points);
        movCloudFull = pointCloud(stateTR.Points);
        
        % randomly downsample
        fixCloud = pcdownsample(fixCloudFull,'random',downPercentage);
        movCloud = pcdownsample(movCloudFull,'random',downPercentage);
        
        %
        % rigid ICP of stateStr to endStr
        %
        rigidRegFname = [resultDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg.mat'];
        
        if ~exist(rigidRegFname,'file')
            disp('Do rigid registration')
            % forward registration
            [tform,movRegCloud,rmse] = pcregistericp(movCloud,fixCloud,'InlierRatio',inlierRatio,'Verbose',true,'MaxIterations',maxIter,'Tolerance',tolVal);
            % inverse registration
            [invTform,invMovRegCloud,invRmse] = pcregistericp(fixCloud,movCloud,'InlierRatio',inlierRatio,'Verbose',true,'MaxIterations',maxIter,'Tolerance',tolVal);
            % initialize forward registration with inverted inverse registration
            [tformB,movRegCloudB,rmseB] = pcregistericp(movCloud,fixCloud,'InlierRatio',inlierRatio,'Verbose',true,'MaxIterations',4*maxIter,'Tolerance',tolVal,'InitialTransform',invert(invTform));
            if rmseB<rmse
                tform=tformB;
                movRegCloud=movRegCloudB;
                rmse=rmseB;
            end
            
            if strcmp(version('-release'),'2022b')
                % MATLAB R2022b: is creating rigidtform3d object, premultiply
                % store also in old format, postmultiplication
                tformPost = rigid3d(tform.R',tform.A(1:3,4)');
                save(rigidRegFname,'tform','tformPost','movRegCloud','rmse');
            else
                % created with old version
                tformPost=tform;
                save(rigidRegFname,'tformPost','movRegCloud','rmse');
            end
        else
            disp('Load rigid registration')
            disp(['load ' rigidRegFname])
            a=load(rigidRegFname);
            movRegCloud=a.movRegCloud;
            rmse=a.rmse;
            
            if isfield(a,'tformPost') 
                tform=a.tformPost;
            else
                tform=a.tform;
                if strcmp(version('-release'),'2022b')
                    % create post transformation
                    tformPre=tform;
                    tformPost = rigid3d(tform.R',tform.A(1:3,4)');
                    save(rigidRegFname,'tform','tformPost','movRegCloud','rmse');
                else 
                    keyboard
                end
            end
            
            redoNR = 0;
            if a.rmse>10
                % relatively high rigid registration error
                % try registration the otherway round to get a better initialization
                [invTform,invMovRegCloud,invRmse] = pcregistericp(fixCloud,movCloud,'InlierRatio',inlierRatio,'Verbose',true,'MaxIterations',maxIter,'Tolerance',tolVal);
                [tform,movRegCloud,rmse] = pcregistericp(movCloud,fixCloud,'InlierRatio',inlierRatio,'Verbose',true,'MaxIterations',4*maxIter,'Tolerance',tolVal,'InitialTransform',invert(invTform));
                if rmse<a.rmse
                    if strcmp(version('-release'),'2022b')
                        % MATLAB R2022b: is creating rigidtform3d object, premultiply
                        % store also in old format, postmultiplication
                        tformPost = rigid3d(tform.R',tform.A(1:3,4)');
                        save(rigidRegFname,'tform','tformPost','movRegCloud','rmse');
                    else
                        % created with old version
                        tformPost=tform;
                        save(rigidRegFname,'tformPost','movRegCloud','rmse');
                    end
                    redoNR = 1;  % redo non-rigid registration
                else
                    keyboard
                end
            end
        end
        
        %%  show rigid result
        disp('Figure   6: rigid registration result')
        figure(6)
        plot3(endTR.Points(:,1),endTR.Points(:,2),endTR.Points(:,3),'b.')
        hold on
        plot3(movRegCloud.Location(:,1),movRegCloud.Location(:,2),movRegCloud.Location(:,3),'r.')
        plot3(stateTR.Points(:,1),stateTR.Points(:,2),stateTR.Points(:,3),'k.')
        hold off
        axis equal
        axis tight
        legend([endStr '-fixed'],[stateStr 'rigid'],[stateStr '-moving'])
        title(['RMSE ' num2str(rmse)])
        view([1 -1 0])
        xlabel('x [mm]')
        ylabel('y [mm]')
        zlabel('z [mm]')
        if ~runThrough
            keyboard
        end
        %% non-rigid ICP registration
        %
        %   Source: structured object with fields -
        %       Source.vertices: V x 3 vertices of template model
        %       Source.faces: F x 3 list of connected vertices.
        %       Source.normals: (Optional) FV x 3 list of surface normals. Make
        %           sure to set Options.normals = 1 if using normals.
        fixMesh.vertices = endTeethTR.Points;
        fixMesh.faces = endTeethTR.ConnectivityList;
        fixMesh.normals = vertexNormal(endTeethTR);
        
        % movMesh: stateTR
        movMesh.vertices = stateTR.Points;
        movMesh.faces = stateTR.ConnectivityList;
        movMesh.normals = vertexNormal(stateTR);
        transfStr = '';
        %
        % forward, endTR is fixed, motion of stateTR
        %
        nrRegFname = [resultDir pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr 'TeethNonRigidReg_' nrRegParStr '.mat'];
        if ~exist(nrRegFname,'file') || redoNR
            % initialize with rigid transformation
            nrOptions.rigidM = tform.T(:,1:3);
            % source is movMesh
            % target is fixMesh
            % ouput are non-rigidly transformed movMesh, transformation X
            [nrTransfPointsFull, X] = nricp(movMesh,fixMesh, nrOptions);
            save(nrRegFname,'nrTransfPointsFull','X');
        else
            disp(['load ' nrRegFname])   % X, nrTransfPointsFull
            load(nrRegFname)
        end
        
        %
        % backward, stateTR is fixed, motion of endTR
        %
        invNrRegFname = [resultDir pat{p}.id '_' partStr '_' endStr 'TeethTo' stateStr 'NonRigidReg_' nrRegParStr '.mat'];
        if ~exist(invNrRegFname,'file') || redoNR
            invM=inv(tform.T);
            nrOptions.rigidM = invM(:,1:3);
            [invNrTransfPointsFull, invX] = nricp(fixMesh,movMesh, nrOptions);
            save(invNrRegFname,'invNrTransfPointsFull','invX');
        else
            disp(['load ' invNrRegFname])  % invX, invNrTransfPointsFull
            load(invNrRegFname)
        end
        
        % rigid transform all points of stateStr to create another triangulation
        rigidTransfPointsFull=transformPointsForward(tform,movCloudFull.Location);
        % create triangulation of transformed moving mesh
        rigidMovingTR = triangulation(stateTR.ConnectivityList,rigidTransfPointsFull);
        % rigid result
        movRegCloud = pointCloud(rigidTransfPointsFull);
        
        if exist(stlFileStateNameTeeth,'file') && endN==refN
            % also rigidly transform teeth of stateTR
            rigidTransfTeethPointsFull=transformPointsForward(tform,stateTeethTR.Points);
            rigidMovingTeethTR = triangulation(stateTeethTR.ConnectivityList,rigidTransfTeethPointsFull);
        end
        
        % backward non-rigid result for motion calculation
        
        % non-rigidly transform all teeth of endTR to stateTR
        invNrMovingTR = triangulation(endTeethTR.ConnectivityList,invNrTransfPointsFull);
        %%
        if ~exist(stlFileStateNameTeeth,'file') && endN==refN
            disp([stlFileStateNameTeeth ' missing'])
            disp(['transfer teeth segmentation and occlusion plane from state ' num2str(endN) ' to ' num2str(n)])
            % transfer teeth segmentation and occPlane
            % from end state N to non-end state
            % based on non-rigid registration
            %
            % extract individual teeth closest after registration
            %
            % show individual tooth segmentation on stateTR
            disp('Figure 901: shows transferred occlusion plane and teeth segmentations')
            figure(901)
            fh=trisurf(stateTR);
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            axis equal tight
            [caz,cel]=view([0 -1 0]);
            lightangle(gca,caz,cel)
            title([stateStr ' ' partStr])
            drawnow
            %
            tmpSTLFname = 'tmp.stl';
            % for occlusion plane
            % find closest point in endTR
            disp('propagate occlusion plane points via non-rigid registration')
            occPlaneEndPoints = [occPlane.P1; occPlane.P2; occPlane.P3];
            for p1=1:3
                [idxV,distV] = knnsearch(endTeethTR.Points,occPlaneEndPoints(p1,:));
                disp(['P' num2str(p1) ' in endTR, max distance ' num2str(max(distV))])
                nrP = invNrMovingTR.Points(idxV,:);   % non-rigid result to stateTR
                [idxV,distV] = knnsearch(stateTR.Points,nrP);
                disp(['P' num2str(p1)  ' in ' stateStr ', max distance ' num2str(max(distV))])
                occPlaneStatePoints(p1,:)=stateTR.Points(idxV,:);
            end
            occPlaneState.P1=occPlaneStatePoints(1,:);
            occPlaneState.P2=occPlaneStatePoints(2,:);
            occPlaneState.P3=occPlaneStatePoints(3,:);
            
            figure(901)
            hold on
            plot3(occPlaneState.P1(1),occPlaneState.P1(2),occPlaneState.P1(3),'ro')
            plot3(occPlaneState.P2(1),occPlaneState.P2(2),occPlaneState.P2(3),'go')
            plot3(occPlaneState.P3(1),occPlaneState.P3(2),occPlaneState.P3(3),'bo')
            occPlaneState.midP = (occPlaneState.P2+occPlaneState.P3)/2;
            plot3(occPlaneState.midP(1),occPlaneState.midP(2),occPlaneState.midP(3),'k*')
            line([occPlaneState.P2(1) occPlaneState.P3(1)],[occPlaneState.P2(2) occPlaneState.P3(2)],[occPlaneState.P2(3) occPlaneState.P3(3)])
            line([occPlaneState.P1(1) occPlaneState.midP(1)],[occPlaneState.P1(2) occPlaneState.midP(2)],[occPlaneState.P1(3) occPlaneState.midP(3)])
            hold off
            
            if ~runThrough
                keyboard
            end
            
            % save in matlab file
            save(occPlaneStateName,'occPlaneState')
            
            % for collecting all teeth
            disp('propagate teeth segmentations via non-rigid registration')
            allPoints = [];
            allConn = [];
            prevNoP = 0;
            for to=1:length(seg{p}.endToothTR)
                % find closest points of given tooth in overall end mesh
                [idxV,distV] = knnsearch(endTeethTR.Points,seg{p}.endToothTR{to}.Points);
                if max(distV)>0
                    disp([num2str(to) ', max distance ' num2str(max(distV)) ' should be zero!'])
                end
                
                % invNrTransfPointsTooth: non-rigid transform endTR to stateTR
                invNrTransfPointsTooth = invNrMovingTR.Points(idxV,:);
                invNrMovingToothTR = triangulation(seg{p}.endToothTR{to}.ConnectivityList,invNrTransfPointsTooth);
                
                % - find nearest points in stateTR
                % find closest points of given tooth in overall state mesh
                % do it with knn (faster), find 2 nearest neighbours
                [idxV,distV] = knnsearch(stateTR.Points, invNrMovingToothTR.Points);
                disp([num2str(to) ', max distance ' num2str(max(distV(:)))])
                
                % find all triangles attached to points
                vertexAttached=vertexAttachments(stateTR,idxV(:));
                vertexList = [];
                for vidx = 1:length(vertexAttached)
                    vertexList = [vertexList vertexAttached{vidx}];
                end
                vertexList = unique(vertexList);
                % create triangulation with all points and
                % selected triangles
                testTR = triangulation(stateTR.ConnectivityList(vertexList,:),stateTR.Points);
                
                % find holes
                % holeCellArray
                [holeCellArray,bounding_triangles,holeLengths] = findTriMeshHoles(testTR.ConnectivityList,testTR.Points);
                [~,maxIdx]=max(holeLengths);
                numPointsOfHole = cellfun(@numel,holeCellArray);
                
                % view holes and collect hole IDs
                holeIDList = [];  % collect ID of hole vertices
                figure(65)
                fh=trisurf(testTR);
                set(fh,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                title('Identify Holes'); hold on; axis equal;
                vertices=testTR.Points;
                for ho = 1:length(holeCellArray)
                    hole = holeCellArray{ho};
                    line(vertices(hole,1),vertices(hole,2),vertices(hole,3),'Color','r')
                    holeIDList = [holeIDList hole'];
                end
                hold off
                holeIDList = unique(holeIDList);
                
                % find all triangles attached to hole points
                vertexHoles=vertexAttachments(stateTR,holeIDList');
                vertexHoleList = [];
                for vidx = 1:length(vertexHoles)
                    vertexHoleList = [vertexHoleList vertexHoles{vidx}];
                end
                vertexAllList = unique([vertexList vertexHoleList]);
                
                % add hole vertices
                testTR = triangulation(stateTR.ConnectivityList(vertexAllList,:),stateTR.Points);
                figure(66)
                fh=trisurf(testTR);
                set(fh,'FaceColor',seg{p}.color{to},'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                title('Filled Holes'); hold on; axis equal;
                hold off
                
                stlwrite(testTR,tmpSTLFname,'text');
                % cleanup by reading it in again
                stateToothTR=stlread(tmpSTLFname);
                
                figure(901)
                hold on
                fh=trisurf(stateToothTR);
                set(fh,'FaceColor',seg{p}.color{to},'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                hold off
                drawnow
                if doUpper
                    outToothFname = [stlFileStateNameToothStr seg{p}.upperToothStrV{to} '.stl'];
                else
                    outToothFname = [stlFileStateNameToothStr seg{p}.lowerToothStrV{to} '.stl'];
                end
                stlwrite(stateToothTR,outToothFname,'text');
                
                % compose from all segmented teeth
                allPoints = [allPoints; stateToothTR.Points];
                % need to add point number offset
                allConn = [allConn; stateToothTR.ConnectivityList+prevNoP];
                prevNoP=prevNoP+size(stateToothTR.Points,1);
            end
            if ~runThrough
                keyboard
            end
            % not existing yet, compose from segmented teeth
            stateTeethTR = triangulation(allConn,allPoints);
            stlwrite(stateTeethTR,stlFileStateNameTeeth,'text');
            
            disp('Figure 902: shows transferred teeth segmentation')
            figure(902)
            fh=trisurf(stateTR);
            set(fh,'FaceColor',BMCblue,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            hold on
            fh=trisurf(stateTeethTR);
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            hold off
            axis equal tight
            [caz,cel]=view([0 -1 0]);
            lightangle(gca,caz,cel)
            legend([stateStr ' ' partStr],'registered teeth','Location','best')
            disp('transfer of teeth and occlusion plane finished')
            break
        end
        
        % forward non-rigid result, stateTR to endTR
        nrMovingTR = triangulation(stateTR.ConnectivityList,nrTransfPointsFull);
        
        %% show registration results
        disp('Figure  40: show side-by-side after rigid registration')
        % show both side by side
        figure(40)
        subplot(1,2,1)
        fh=trisurf(endTR);
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        axis tight
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title([pat{p}.id ' ' partStr ' ' endStr])
        [caz,cel]=view([1 -1 0]);
        lightangle(gca,caz,cel)
        material(materialV)
        
        % set limits for rendering/movies where lower jaw is rotated
        if doUpper
            xlimRot=xlimSave;
            ylimRot=ylimSave;
            zlimRot=zlimSave;
            xlimTeethRot=xlimSave;
            ylimTeethRot=ylimSave;
            zlimTeethRot=zlimSave;
        else
            % turn x and z round as lower is rotated
            xlimRot=xlimSave;
            ylimRot=[-ylimSave(2) -ylimSave(1)];
            zlimRot=zlimSave;
            %
            ylimTeethRot=[-ylimSave(2) -ylimSave(1)];
            zlimTeethRot=zlimSave;
            xlimTeethRot=xlimSave;
        end
        xlim(xlimSave)
        ylim(ylimSave)
        zlim(zlimSave)
        %
        subplot(1,2,2)
        fh=trisurf(rigidMovingTR);
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        axis tight
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title([pat{p}.id ' ' partStr ' rigid ' stateStr])
        [caz,cel]=view([1 -1 0]);
        lightangle(gca,caz,cel)
        material(materialV)
        
        % make it same limits
        xlim(xlimSave)
        ylim(ylimSave)
        zlim(zlimSave)
        if saveFigs
            print(40,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_sidebyside.png'])
            subplot(1,2,1)
            view([0 -1 0]); 
            lightangle(gca,-caz,cel)
            subplot(1,2,2)
            view([0 -1 0]);
            lightangle(gca,-caz,cel)
            print(40,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_sidebysideY.png'])
        end
        
        disp('Figure  50: show overlay after rigid registration')
        % show overlayed
        figure(50)
        fh=trisurf(endTR);
        set(fh,'FaceColor',BMCblue,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        % stick to fixed image axis
        xlim(xlimSave)
        ylim(ylimSave)
        zlim(zlimSave)
        
        xlabel('x')
        ylabel('y')
        zlabel('z')
        [caz,cel]=view([1 -1 0]);
        lightangle(gca,caz,cel)
        material(materialV)
        %
        hold on
        fh=trisurf(rigidMovingTR);
        if doEndSimReg
            set(fh,'FaceColor',BMCred,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            title([pat{p}.id ' ' partStr ', ' endStr ' (blue) rigid ' stateStr ' (red)'])
        else
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            title([pat{p}.id ' ' partStr ', ' endStr ' (blue) rigid ' stateStr ' (white)'])
        end
        hold off
        if saveFigs
            print(50,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_both.png'])
            % add second light source
            lightangle(gca,-1*caz,cel)
            % from top
            view([0 -1 0])
            print(50,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_both_view3.png'])
            % from front
            view([0 0 1])
            print(50,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_both_view2.png'])
        end
        if ~runThrough
            keyboard
        end
        %% save individual images to later make a movie
        
        if (saveFigs || saveMainFigs) && n==refN-1 
            % save also end state, with out difference teeth and gum color
            figure(51)
            fh=trisurf(endTR);
            if doUpper==0
                rotate(fh,[1 0 0],180);
            end
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',1,'FaceLighting',fLighting)
            axis equal
            axis tight
            % from top
            if doUpper
                [caz,cel]=view([0 -1 0]);
            else
                [caz,cel]=view([0 1 0]);
            end
            lightangle(gca,caz,cel)
            grid off
            axis off
            % set background color, lighter BMCblue
            set(gcf,'color',lightBMCblue)
            
            xlim(xlimRot)
            ylim(ylimRot)
            zlim(zlimRot)
            
            material(materialV)
            set(51, 'InvertHardCopy', 'off');
            print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) 'XTo' endStr '_rigidReg.png'])
        end
        
        figure(51)
        fh=trisurf(endTR);
        if doUpper==0
            rotate(fh,[1 0 0],180);
        end
        set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',1,'FaceLighting',fLighting)
        hold on
        fh=trisurf(endTeethTR);
        if doUpper==0
            rotate(fh,[1 0 0],180);
        end
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        hold off
        axis equal
        axis tight
        % from top
        if doUpper
            [caz,cel]=view([0 -1 0]);
        else
            [caz,cel]=view([0 1 0]);
        end
        lightangle(gca,caz,cel)
        grid off
        axis off
        % set background color, lighter BMCblue
        set(gcf,'color',lightBMCblue)
        
        xlim(xlimRot)
        ylim(ylimRot)
        zlim(zlimRot)
        
        material(materialV)
        if (saveFigs || saveMainFigs) && n==refN-1 
            % save also end state
            set(51, 'InvertHardCopy', 'off');
            print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) 'XTo' endStr '_rigidReg_gum.png'])
        end
        if saveFigs
            % output fixed image
            % keep background
            set(51, 'InvertHardCopy', 'off');
            % end is always in the middle of movie
            
            print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_002.png'])
            if makeClean
                hold on
                if doUpper==0
                    hold off
                    % no turning
                    figure(51)
                    fh=trisurf(endTR);
                    set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',1,'FaceLighting',fLighting)
                    hold on
                    fh=trisurf(endTeethTR);
                    set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                    hold off
                    axis equal
                    axis tight
                    [caz,cel]=view([0 -1 0]);
                    lightangle(gca,caz,cel)
                    grid off
                    axis off
                    % set background color, lighter BMCblue
                    xlim(xlimSave)
                    ylim(ylimSave)
                    zlim(zlimSave)
                    
                    material(materialV)
                    hold on
                    plot3(occPlane.P1(1),occPlane.P1(2),occPlane.P1(3),'r*','MarkerSize',ms,'LineWidth',lw)
                    plot3(occPlane.P2(1),occPlane.P2(2),occPlane.P2(3),'g*','MarkerSize',ms,'LineWidth',lw)
                    plot3(occPlane.P3(1),occPlane.P3(2),occPlane.P3(3),'b*','MarkerSize',ms,'LineWidth',lw)
                    occPlane.midP = (occPlane.P2+occPlane.P3)/2;
                    line([occPlane.P2(1) occPlane.P3(1)],[occPlane.P2(2) occPlane.P3(2)],[occPlane.P2(3) occPlane.P3(3)],'LineWidth',lw,'color','k')
                    line([occPlane.P1(1) occPlane.midP(1)],[occPlane.P1(2) occPlane.midP(2)],[occPlane.P1(3) occPlane.midP(3)],'LineWidth',lw,'color','k')
                else
                    plot3(occPlane.P1(1),occPlane.P1(2),occPlane.P1(3),'r*','MarkerSize',ms,'LineWidth',lw)
                    plot3(occPlane.P2(1),occPlane.P2(2),occPlane.P2(3),'g*','MarkerSize',ms,'LineWidth',lw)
                    plot3(occPlane.P3(1),occPlane.P3(2),occPlane.P3(3),'b*','MarkerSize',ms,'LineWidth',lw)
                    occPlane.midP = (occPlane.P2+occPlane.P3)/2;
                    line([occPlane.P2(1) occPlane.P3(1)],[occPlane.P2(2) occPlane.P3(2)],[occPlane.P2(3) occPlane.P3(3)],'LineWidth',lw,'color','k')
                    line([occPlane.P1(1) occPlane.midP(1)],[occPlane.P1(2) occPlane.midP(2)],[occPlane.P1(3) occPlane.midP(3)],'LineWidth',lw,'color','k')
                end
                hold off
                set(51, 'InvertHardCopy', 'on');
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_occPlane_gum.png'])
           
            else
                if doUpper
                    [caz,cel]=view([0 0 1]);
                    % add second light
                    %lightangle(gca,-1*caz,30)
                    lightangle(gca,caz,cel)
                else
                    [caz,cel]=view([0 0 -1]);
                    lightangle(gca,caz,cel)
                end
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_view2_002.png'])
            end
        end
             
        if  ((saveFigs || saveMainFigs) || redoMovie ) && n==refN-1
            % show teeth rendering for movie
            figure(51)
            fh=trisurf(endTeethTR);
            if doUpper==0
                rotate(fh,[1 0 0],180);
            end
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            axis equal
            axis tight
            % from top
            if doUpper
                [caz,cel]=view([0 -1 0]);
            else
                [caz,cel]=view([0 1 0]);
            end
            lightangle(gca,caz,cel)
            grid off
            axis off
            xlim(xlimTeethRot)
            ylim(ylimTeethRot)
            zlim(zlimTeethRot)
            if redoMovie
                % set background color, light BMCblue
                set(gcf,'color',gumColor)
            else
                % set background color, light BMCblue
                set(gcf,'color',lightBMCblue)
            end
            material(materialV)
            set(51, 'InvertHardCopy', 'off');
            
            if redoMovie
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) 'XTo' endStr '_teeth_rigidReg_gumColor.png'])
            else
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) 'XTo' endStr '_teeth_rigidReg.png'])
            end
        end
        
        if  n==1 && endN==refN
            % show rendering of achieved, planned, and first teeth
            figure(520)
            if ~doEndSimReg
                fh=trisurf(endTeethTR);
                % final measured endTR is BMCblue
                set(fh,'FaceColor',BMCblue,'EdgeColor','none','FaceAlpha',1,'FaceLighting',fLighting)
            end
            hold on
            % registered teet from state n or simState n
            fh=trisurf(rigidMovingTeethTR);
            if doEndSimReg
                % rigid aligned simState n is BMCred
                set(fh,'FaceColor',BMCred,'EdgeColor','none','FaceAlpha',1,'FaceLighting',fLighting)
            else
                % rigid aligned state n is BMCorange
                set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',1,'FaceLighting',fLighting)
            end
            axis equal
            axis tight
        
            % from top
            [caz,cel]=view([0 -1 0]);
            if ~doEndSimReg
                lightangle(gca,caz,cel)
            end
            grid off
            axis off
            xlim(xlimSave)
            ylim(ylimSave)
            zlim(zlimSave)
            material(materialV)
            if saveFigs
                print(520,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg_All.png'])
            end
            
            % add a length bar on upper figures
            if doUpper && doEndSimReg
                % only on final
                minValues=min(endTeethTR.Points);
                maxValues=max(endTeethTR.Points);
                meanValues=mean(endTeethTR.Points);
                centerValues = minValues+maxValues;
                hold on
                line([-5 5]+meanValues(1),[0 0]+centerValues(2),[8 8]+minValues(3),'LineWidth',5,'Color','k')
                text(meanValues(1),centerValues(2),4+minValues(3),'10 mm','HorizontalAlignment','center','Color','k','FontSize',3*fs)
                if saveFigs
                    print(520,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg_All_SB.png'])
                end
            end
            
            if saveGAbstractFigs && ~doEndSimReg
                % in Week 0 in blue, End in red. No Simulation
                figure(521), clf
                fh=trisurf(endTeethTR);
                % final measured endTR is BMCred
                set(fh,'FaceColor',BMCred,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                hold on
                % registered teet from state n
                fh=trisurf(rigidMovingTeethTR);
                % rigid aligned state n is BMCblue
                set(fh,'FaceColor',BMCblue,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                axis equal
                axis tight
                % from top
                [caz,cel]=view([0 -1 0]);
                lightangle(gca,caz,cel)
                grid off
                axis off
                xlim(xlimSave)
                ylim(ylimSave)
                zlim(zlimSave)
                material(materialV)
                print(521,'-dpng','tmpFlipX.png','-r300')
                if doUpper
                    toXFlipI=imread('tmpFlipX.png');
                    imwrite(toXFlipI,[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg_All_redblue.png'])
                else
                    % needs rotating and flipping, as rotating fh makes
                    % overlay bad
                    toXFlipI=imrotate(imread('tmpFlipX.png'),180);
                    imwrite(fliplr(toXFlipI),[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg_All_redblue.png'])
                end
                        
            end
        end
        
        % show gum, all in white
        disp('Figure  51: rigidly aligned state for movie')
        figure(51),clf
        fh=trisurf(rigidMovingTR);
        if doUpper==0
            rotate(fh,[1 0 0],180);
        end
        set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        % from top
        if doUpper
            [caz,cel]=view([0 -1 0]);
        else
            [caz,cel]=view([0 1 0]);
        end
        lightangle(gca,caz,cel)
        if redoMovie
            % set background color, light BMCblue
            set(gcf,'color',gumColor)
        else
            set(gcf,'color',lightBMCblue)
        end
        set(51, 'InvertHardCopy', 'off');
        material(materialV)
        axis off
        
        % make it same limits
        xlim(xlimRot)
        ylim(ylimRot)
        zlim(zlimRot)
        
        if endN==refN 
            % save State n registered to EndN state for making moving over all end states
            if redoMovie
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_rigidReg_gumColor.png'])
            elseif saveFigs
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_rigidReg.png'])
            end
            % also show teeth
            if exist(stlFileStateNameTeeth,'file')
                figure(51),clf
                fh=trisurf(rigidMovingTeethTR);
                if doUpper==0
                    rotate(fh,[1 0 0],180);
                end
                set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                axis equal
                axis off
                % from top
                if doUpper
                    [caz,cel]=view([0 -1 0]);
                else
                    [caz,cel]=view([0 1 0]);
                end
                lightangle(gca,caz,cel)
                if redoMovie
                    % set background color to gum color
                    set(gcf,'color',gumColor)
                else
                    % set background color, light BMCblue
                    set(gcf,'color',lightBMCblue)
                end
                set(51, 'InvertHardCopy', 'off');
                material(materialV)
                
                xlim(xlimTeethRot)
                ylim(ylimTeethRot)
                zlim(zlimTeethRot)
                if redoMovie
                    print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg_gumColor.png'])
                elseif saveFigs
                    print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg.png'])
                end
                if redoMovie && n==refN-1
                    % to make movie in unix
                    gifCmd =['convert -adjoin -delay 100 -loop 0 ' regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) '?To' endStr '_teeth_rigidReg.png ' regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) 'sTo' endStr '_teeth_rigidReg_movie.gif'];
                    if isunix
                        unix(gifCmd)
                    else
                        disp(gifCmd)
                    end
                end
            end
            if redoMovie && n==refN-1
                % to make movie in unix
                gifCmd =['convert -adjoin -delay 100 -loop 0 ' regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) '?To' endStr '_rigidReg.png ' regFigsDir pat{p}.id '_' partStr '_' stateStr(1:end-1) 'sTo' endStr '_rigidReg_movie.gif'];
                if isunix
                    unix(gifCmd)
                else
                    disp(gifCmd)
                end
            end
        end
        if saveFigs
            % state is 001, simulationEnd is 003
            if doEndSimReg
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_003.png'])
                if doUpper
                    [caz,cel]=view([0 0 1]);
                    % add second light
                    lightangle(gca,caz,cel)
                else
                    [caz,cel]=view([0 0 -1]);
                    lightangle(gca,caz,cel)
                end
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_view2_003.png'])
            else
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_001.png'])
                if doUpper
                    [caz,cel]=view([0 0 1]);
                    % add second light
                    lightangle(gca,caz,cel)
                else
                    [caz,cel]=view([0 0 -1]);
                    lightangle(gca,caz,cel)
                end
                print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_view2_001.png'])
            end
            % to make movie in unix
            disp(['convert -adjoin -delay 100 -loop 0 ' regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_0??.png ' regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_movie.gif'])
            disp(['convert -adjoin -delay 100 -loop 0 ' regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_view2_0??.png ' regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_view2_movie.gif'])
        end
        
        if ~runThrough
            keyboard
        end
        %% show distance to closest point on reference mesh
        
        % determine closest point from endTR to movRegCloud (rigidly aligned stateTR) 
        % and Euclidean distance to it
        disp('determine distance to closest point after rigid registration')
        % knnsearch much faster than pdist2
        [closestIdx,EuclD] = knnsearch(movRegCloud.Location,endTR.Points);
                            
        % determine component errors
        % endTR - rigid(movTR)
        componentD = endTR.Points-movRegCloud.Location(closestIdx,:);
        
        %% show distance between rigidly registered surfaces
        errorStr = '';
        
        for i=1:4
            if i==4
                % get also closestIdx
                minDistancesFix2Mov = EuclD;
                disp('Figure   7: Euclidean distance to closest point')
            else
                if doEndSimReg
                    % from endTR to endSimTR
                    minDistancesFix2Mov = -componentD(:,i)';
                    dispStr = [endStr '->' stateStr];
                else
                    % from stateTR to endTR
                    minDistancesFix2Mov = componentD(:,i)';
                    dispStr = [stateStr '->' endStr];
                end
            end
            figure(7)
            fh=trisurf(endTR);
            % set error as color
            set(fh,'CData',1000*minDistancesFix2Mov,'FaceColor','interp','EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            axis equal tight
            colorbar
            xlabel('x [mm]')
            ylabel('y [mm]')
            zlabel('z [mm]')
            [caz,cel]=view([1 -1 0]);
            lightangle(gca,caz,cel)
            material(materialV)
            
            meanD = mean(minDistancesFix2Mov);
            stdD = std(minDistancesFix2Mov);
            errorStr = [errorStr num2str(1000*meanD,'%.0f') '$\pm$' num2str(1000*stdD,'%.0f') ' & '];
            if i==4
                % show only up to inlier ratio distance
                caxis(1000*[0 1])
                title([stateStr 'To' endStr ' EuclideanD, mean ' num2str(1000*meanD,'%.0f') ' um'])
                colormap('parula')
            else
                caxis(1000*[-1 1])
                title([dispStr ' in ' coordStrV{i} ', mean ' num2str(1000*meanD,'%.0f') '+-' num2str(1000*stdD,'%.0f') ' um'])
                colormap(redblue)
            end
            if saveFigs
                if i==4
                    print(7,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_EuclDist' num2str(i) '.png'])
                    
                    maxErr=prctile(minDistancesFix2Mov,100*inlierRatio);
                    caxis([0 1000*maxErr])
                    print(7,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_EuclDist' num2str(i) 'Inliers.png'])
                    
                else
                    print(7,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr 'RigidReg_DistX' num2str(i) '.png'])
                end
            end
        end
        % list results
        disp(' & & & $x$ & $y$ & $z$ & 3D & RMSE \\')
        disp([pat{p}.id ' & ' extraOutStr ' ' partStr ' & ' stateStr 'To' endStr ' & ' errorStr num2str(1000*rmse,'%.0f') ' um \\'])
        if ~runThrough
            keyboard
        end
        %% non-rigid registration, all aligning, more important is motion
        
        disp('Figure 500: overlay of non-rigid result')
        figure(500)
        fh=trisurf(endTeethTR);
        teethStr = 'Teeth';
        set(fh,'FaceColor',BMCblue,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
        axis equal
        % stick to fixed image axis
        xlim(xlimSave)
        ylim(ylimSave)
        zlim(zlimSave)
        
        xlabel('x')
        ylabel('y')
        zlabel('z')
        [caz,cel]=view([1 -1 0]);
        lightangle(gca,caz,cel)
        material(materialV)
        %
        hold on
        fh=trisurf(nrMovingTR);
        if doEndSimReg
            set(fh,'FaceColor',BMCred,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            title([pat{p}.id ' ' partStr ', ' endStr teethStr ' (blue) ' transfStr stateStr ' (red) non-rigid'])
        else
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            title([pat{p}.id ' ' partStr ', ' endStr teethStr ' (blue) ' transfStr stateStr ' (white) non-rigid'])
        end
        hold off
        if saveFigs
            print(500,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidReg_both.png'])
            view([0 0 1])
            lightangle(gca,-1*caz,cel)
            print(500,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidReg_both_view2.png'])
        end
        
        %% motion from non-rigid registration
        % rigidly transform fixed endTR to moving
        invRigidTransfPointsFull=transformPointsInverse(tform,endTeethTR.Points);
        % rigidly transformed End to stateTR - non-rigid transformed
        % End to stateTR
        if doEndSimReg
            % from rigid(endTR) to endSimTR
            nrMovDisp = invNrTransfPointsFull - invRigidTransfPointsFull;
            dispStr = [endStr '->' stateStr];
            tableStr = [endStr '$\rightarrow$' stateStr];
        else
            % from stateTR to rigid(endTR)
            nrMovDisp = invRigidTransfPointsFull - invNrTransfPointsFull;
            dispStr = [stateStr '->' endStr];
            tableStr = [stateStr '$\rightarrow$' endStr];
        end
        errorStr = '';
        % i=1-3: motion compontents in direction 1-3
        % i=4:5: Euclidean distance
        for i=[2 5] %1:5
            if i<4
                disp(['Figure ' num2str(400+i) ': motion component in direction ' num2str(i)])
            else
                disp(['Figure ' num2str(400+i) ': motion Euclidean distance'])
            end
            figure(400+i)
            % nonrigidly transformed
            fh=trisurf(endTeethTR);
            if i>3
                nrDist=sqrt(sum(nrMovDisp.^2,2));
            else
                nrDist=nrMovDisp(:,i);
            end
            
            set(fh,'CData',1000*nrDist,'FaceColor','interp','EdgeColor','none','EdgeAlpha',fAlpha,'FaceAlpha',fAlpha,'FaceLighting',fLighting)
            axis equal tight
            colorbar
            if i==5 || i==2
                [caz,cel]=view([0 -1 0]);
                axis tight
            else
                xlabel('x [mm]')
                ylabel('y [mm]')
                zlabel('z [mm]')
                [caz,cel]=view([1 -1 0]);
            end
            % same size
            xlim(occPlane.midP(1)+xlimMidP)
            ylim(occPlane.midP(2)+ylimMidP)
            zlim(occPlane.midP(3)+zlimMidP)
            lightangle(gca,caz,cel)
            material(materialV)
            
            meanD = mean(nrDist);
            stdD = std(nrDist);
            errorStr = [errorStr num2str(1000*meanD,'%.0f') '$\pm$' num2str(1000*stdD,'%.0f') ' & '];
            if i>3
                % show only up to inlier ratio distance
                caxis(1000*[0 1])
                if i==4
                    title([stateStr 'To' endStr ' EuclideanD, mean ' num2str(1000*meanD,'%.0f') ' um'])
                end
                %colStr = '';     % parula
                %colStr = '_BCYR';
                %colStr = '_intBCYR';
                %colStr = '_blindBCYR';
                %colStr = '_BCY';
                %colStrV = {'_BCYR','_intBCYR','_blindBCYR','_BCY','_intBCY','_blindBCY'};
                if doEndSimReg
                    if saveGAbstractFigs
                        colStrV = {'_intYCB'};
                    else
                        colStrV = {'_intBCYR'};
                    end
                else
                    colStrV = {'_intBCY'};
                end
                % create colormaps
                for c=1:length(colStrV)
                    colStr = colStrV{c};
                    if strcmp(colStr,'_BCYR')
                        % make it a uniform colorbar
                        colM(1,:)=[0 0 1]; % blue
                        colM(2,:)=[0 1 1]; % cyan
                        colM(3,:)=[1 1 0]; % yellow
                        colM(4,:)=[1 0 0]; % red
                        colP=[1 85 170 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_intBCYR')
                        % make it a uniform colorbar
                        % increase illumination
                        colM(1,:)=[0 0 0.25]; % blue
                        colM(2,:)=[0 0.5 0.5]; % cyan
                        colM(3,:)=[0.75 0.75 0]; % yellow
                        colM(4,:)=[1 0 0]; % red
                        colP=[1 85 170 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_intYCB')
                        % make it a uniform colorbar
                        colM(1,:)=[0.75 0.75 0]; % yellow
                        colM(2,:)=[0 0.5 0.5]; % cyan
                        colM(3,:)=[0 0 0.25]; % blue
                        colP=[1 128 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_intCB')
                        % make it a uniform colorbar
                        colM(1,:)=[0 0.5 0.5]; % cyan
                        colM(2,:)=[0 0 0.25]; % blue
                        colP=[1 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_blindBCYR')
                        % colorblind friendly
                        colM(1,:)=[78 121 197]/255; % blue
                        colM(2,:)=[105 177 144]/255; % cyan
                        colM(3,:)=[221 170 51]/255; % yellow
                        colM(4,:)=[218 34 34]/255; % red
                        colP=[1 85 170 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_intBCY')
                        % BCY with increasing illumination
                        colM(1,:)=[0 0 0.33]; % blue
                        colM(2,:)=[0 0.67 0.67]; % cyan
                        colM(3,:)=[1 1 0]; % yellow
                        colP=[1 128 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_BCY')
                        % make it a uniform colorbar
                        colM(1,:)=[0 0 1]; % blue
                        colM(2,:)=[0 1 1]; % cyan
                        colM(3,:)=[1 1 0]; % yellow
                        colP=[1 128 256];
                        colorMV{c}=ownColors(colM,colP);
                    elseif strcmp(colStr,'_blindBCY')
                        % colorblind friendly
                        colM(1,:)=[78 121 197]/255; % blue
                        colM(2,:)=[105 177 144]/255; % cyan
                        colM(3,:)=[221 170 51]/255; % yellow
                        colP=[1 128 256];
                        colorMV{c}=ownColors(colM,colP);
                    else
                        colorMV{c} = colormap('parula');
                    end
                end
            else
                caxis(400*[-1 1])
                if i==1
                    title([dispStr ' in ' coordStrV{i} ', mean ' num2str(1000*meanD,'%.0f') '+-' num2str(1000*stdD,'%.0f') ' um'])
                end
                if i==2
                    % make it a uniform colorbar
                    clear colM
                    % Basel mint, very weak
                    colM(1,:)=[165 215 210]/255; % mint
                    colM(2,:)=[1 1 1]; % white
                    colM(3,:)=[242 181 153]/255; % complement
                    % stronger mint, increase saturation
                    hsv1=rgb2hsv(colM(1,:));
                    hsv1(2)=1.0;   % increase saturation
                    hsv1(3)=hsv1(3)*0.7;    % decrease value to darken
                    colM(1,:)=hsv2rgb(hsv1);
                    hsv1=rgb2hsv(colM(3,:));
                    hsv1(2)=1.0;
                    hsv1(3)=hsv1(3)*0.7;  % decrease value to make darker
                    colM(3,:)=hsv2rgb(hsv1);
                    % make color map
                    colP=[1 128 256];
                    colorMmint=ownColors(colM,colP);
                    colormap(colorMmint)
                else
                    colormap(redblue)
                end
            end
            if saveFigs || (saveMainFigs && (i==5 || i==2)) || (saveGAbstractFigs && i==5)
                if i==4
                    for c=1:length(colStrV)
                        colStr = colStrV{c};
                        figure(400+i)
                        colormap(colorMV{c})
                        pause(1)
                        print(400+i,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist' num2str(i) colStr '.png'])
                    end
                elseif i==5
                    axis off
                    grid off
                    colorbar off
                    for c=1:length(colStrV)
                        colStr = colStrV{c};
                        figure(400+i)
                        colormap(colorMV{c})
                        pause(1)
                        % Fig 2 in manuscript
                        print(400+i,'-dpng','tmpFlipX.png')
                        toXFlipI=imread('tmpFlipX.png');
                        if doUpper==0 && doXFlipLowerDisp
                            % flip in x and rotate 180, as flip in y
                            if doEndSimReg
                                imwrite(flipud(toXFlipI),[regFigsDir 'Fig7_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist_clean' colStr '.png'])
                            else
                                imwrite(flipud(toXFlipI),[regFigsDir 'Fig2_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist_clean' colStr '.png'])
                            end
                        else
                            if doEndSimReg
                                imwrite(toXFlipI,[regFigsDir 'Fig7_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist_clean' colStr '.png'])
                            else
                                imwrite(toXFlipI,[regFigsDir 'Fig2_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist_clean' colStr '.png'])
                            end
                        end
                        
                        % make plain colorbar ratio 18x1
                        clear colorBarI
                        for i1=1:3
                            colorBarI(:,:,i1)=repmat(colorMV{c}(:,i1),[1 round(256/18)])';
                        end
                        if doEndSimReg
                            imwrite(colorBarI,[regFigsDir 'Fig7_colorbar' colStr '.png'])
                        else
                            imwrite(colorBarI,[regFigsDir 'Fig2_colorbar' colStr '.png'])
                        end
                        
                    end
                    if doUpper && (endN==refN || endN==5)
                        for c=1:length(colStrV)
                            % add a length bar
                            minValues=min(endTeethTR.Points);
                            maxValues=max(endTeethTR.Points);
                            meanValues=mean(endTeethTR.Points);
                            centerValues = minValues+maxValues;
                            hold on
                            line([-5 5]+meanValues(1),[0 0]+centerValues(2),[8 8]+minValues(3),'LineWidth',5,'Color','k')
                            text(meanValues(1),centerValues(2),4+minValues(3),'10 mm','HorizontalAlignment','center','Color','k','FontSize',3*fs)
                            hold off
                            pause(1)
                            print(400+i,'-dpng','tmpFlipX.png')
                            toXFlipI=imread('tmpFlipX.png');
                            if doEndSimReg
                                imwrite(toXFlipI,[regFigsDir 'Fig7_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist_clean' colStr '_LB.png'])
                            else
                                imwrite(toXFlipI,[regFigsDir 'Fig2_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_EuclDist_clean' colStr '_LB.png'])
                            end
                        end
                    end
                elseif i==2
                    axis off
                    grid off
                    print(400+i,'-dpng','tmpFlipX.png')
                    toXFlipI=imread('tmpFlipX.png');
                    if doUpper==0 && doXFlipLowerDisp
                        % flip in x and rotate 180, as flip in y
                        imwrite(flipud(toXFlipI),[regFigsDir 'Fig3_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_DistX' num2str(i) '_mint.png'])
                    else
                        imwrite(toXFlipI,[regFigsDir 'Fig3_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_DistX' num2str(i) '_mint.png'])
                    end
                    colorbar off
                    print(400+i,'-dpng','tmpFlipX.png')
                    toXFlipI=imread('tmpFlipX.png');
                    if doUpper==0 && doXFlipLowerDisp
                        % flip in x and rotate 180, as flip in y
                        imwrite(flipud(toXFlipI),[regFigsDir 'Fig3_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_DistX' num2str(i) '_clean_mint.png'])
                    else
                        imwrite(toXFlipI,[regFigsDir 'Fig3_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_DistX' num2str(i) '_clean_mint.png'])
                    end
                    if doUpper && (endN==refN || endN==5)
                        % add a length bar
                        minValues=min(endTeethTR.Points);
                        maxValues=max(endTeethTR.Points);
                        meanValues=mean(endTeethTR.Points);
                        centerValues = minValues+maxValues;
                        hold on
                        line([-5 5]+meanValues(1),[0 0]+centerValues(2),[8 8]+minValues(3),'LineWidth',5,'Color','k')
                        text(meanValues(1),centerValues(2),4+minValues(3),'10 mm','HorizontalAlignment','center','Color','k','FontSize',3*fs)
                        hold off
                        pause(1)
                        print(400+i,'-dpng','tmpFlipX.png')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(toXFlipI,[regFigsDir 'Fig3_' pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_DistX' num2str(i) '_clean_mint_LB.png'])
                    end
                    % make plain colorbar ratio 18x1
                    clear colorBarI
                    for i1=1:3
                        colorBarI(:,:,i1)=repmat(colorMmint(:,i1),[1 round(256/18)])';
                    end
                    imwrite(colorBarI,[regFigsDir 'Fig3_colorbar_mint.png'])
                else
                    if saveFigs
                        print(400+i,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' transfStr stateStr 'To' endStr teethStr 'NonRigidRegInv_DistX' num2str(i) '.png'])
                    end
                end
            end
        end
        disp([pat{p}.id ' & ' extraOutStr ' ' partStr ' & ' transfStr dispStr teethStr ' NR & ' errorStr  ' \\'])
        %
        % save the motion between rigid and non-rigid alignment
        regSave{p}.meanD{doEndSimReg+1} = nrMovDisp;
        
        % show relative achievement
        if doEndSimReg==1 && length(doEndSimRegV)>1
            % show vector on endTR
            nrMovDispT1 = regSave{p}.meanD{1};  % from stateTR to rigid(endTR)
            nrMovDispT3 = nrMovDispT1 + regSave{p}.meanD{2}; % from stateTR to rigid(endTR) to sim N TR
            
            meanT3 = mean(nrMovDispT3);
            stdT3 = std(nrMovDispT3);
            errT3 = sqrt(sum(nrMovDispT3.^2,2));
            errorStrT3 = '';
            for i=1:3
                errorStrT3 = [errorStrT3 num2str(1000*meanT3(i),'%.0f') '$\pm$' num2str(1000*stdT3(i),'%.0f') ' & '];
            end
            errorStrT3 = [errorStrT3 num2str(1000*mean(errT3),'%.0f') '$\pm$' num2str(1000*std(errT3),'%.0f') ' & '];
            
            disp([pat{p}.id ' & ' extraOutStr ' ' partStr ' & ' transfStr 'Start2Sim' endStr teethStr ' NR & ' errorStrT3  ' \\'])
        end
        if ~runThrough
            keyboard
        end
        
        
        
        %% evaluation per individual tooth
        if ~isempty(seg{p}.toothStrV) 
            
            noT = length(seg{p}.endToothTR);  % number of teeth
            noS1=round(sqrt(noT));   % number of subplots
            noS2=ceil(noT/noS1);
            
            % determine tooth coordinate system
            % teeth are entered in order
            % labels provide anatomical grouping (1-4) and order (incisor 1 - molar 8)
            legStr = {};
            if noT>0
                toothCenterV=zeros(3,noT);  % tooth center location
                toothGroupV=zeros(1,noT);   % tooth group label (1-4)
                toothNumberV=zeros(1,noT);  % tooth order (incisor 1 - molar 8)
                for to=1:noT
                    toothCenterV(:,to)=mean(seg{p}.endToothTR{to}.Points);
                    toothGroupV(to) = str2num(seg{p}.toothStrV{to}(1));
                    toothNumberV(to) = str2num(seg{p}.toothStrV{to}(2));
                end
                % show center locations, labels and teeth where occlusion
                % plane was defined
                disp('Figure 329: local coordinate system of crowns')
                figure(329)
                plot3(toothCenterV(1,:),toothCenterV(2,:),toothCenterV(3,:),'rx-')
                legStr = [legStr 'teethCenters'];
                hold on
                for to=1:noT
                    text(toothCenterV(1,to),toothCenterV(2,to),toothCenterV(3,to),seg{p}.toothStrV{to},'FontSize',fs)
                    % show teeth in neighbourhood of occlusion plane points
                    if doUpper
                        if strcmp(seg{p}.toothStrV{to},'11') || strcmp(seg{p}.toothStrV{to},'17')
                            plot3(seg{p}.endToothTR{to}.Points(:,1),seg{p}.endToothTR{to}.Points(:,2),seg{p}.endToothTR{to}.Points(:,3),'c.')
                            legStr = [legStr seg{p}.toothStrV{to}];
                        end
                        if strcmp(seg{p}.toothStrV{to},'21') || strcmp(seg{p}.toothStrV{to},'27')
                            plot3(seg{p}.endToothTR{to}.Points(:,1),seg{p}.endToothTR{to}.Points(:,2),seg{p}.endToothTR{to}.Points(:,3),'y.')
                            legStr = [legStr seg{p}.toothStrV{to}];
                        end
                    else
                        if strcmp(seg{p}.toothStrV{to},'41') || strcmp(seg{p}.toothStrV{to},'47')
                            plot3(seg{p}.endToothTR{to}.Points(:,1),seg{p}.endToothTR{to}.Points(:,2),seg{p}.endToothTR{to}.Points(:,3),'c.')
                            legStr = [legStr seg{p}.toothStrV{to}];
                        end
                        if strcmp(seg{p}.toothStrV{to},'31') || strcmp(seg{p}.toothStrV{to},'37')
                            plot3(seg{p}.endToothTR{to}.Points(:,1),seg{p}.endToothTR{to}.Points(:,2),seg{p}.endToothTR{to}.Points(:,3),'y.')
                            legStr = [legStr seg{p}.toothStrV{to}];
                        end
                    end
                end
                hold off
                axis equal
                xlabel('x')
                ylabel('y')
                zlabel('z')
                grid on
                view([1 -1 0]);
                title([endStr ' ' partStr])
                % show occlusion plane
                if ~isempty(occPlane.P1)
                    hold on
                    plot3(occPlane.P1(1),occPlane.P1(2),occPlane.P1(3),'ro')
                    plot3(occPlane.P2(1),occPlane.P2(2),occPlane.P2(3),'go')
                    plot3(occPlane.P3(1),occPlane.P3(2),occPlane.P3(3),'bo')
                    occPlane.midP = (occPlane.P2+occPlane.P3)/2;
                    plot3(occPlane.midP(1),occPlane.midP(2),occPlane.midP(3),'k*')
                    line([occPlane.P2(1) occPlane.P3(1)],[occPlane.P2(2) occPlane.P3(2)],[occPlane.P2(3) occPlane.P3(3)])
                    line([occPlane.P1(1) occPlane.midP(1)],[occPlane.P1(2) occPlane.midP(2)],[occPlane.P1(3) occPlane.midP(3)])
                    hold off
                    if doUpper
                        legStr = [legStr 'P1','P2 (17)','P3 (27)','midP'];
                    else
                        legStr = [legStr 'P1','P2 (47)','P3 (37)','midP'];
                    end
                end
                % define tooth coordinate system
                toothVectorV = zeros(3,noT);
                if doUpper
                    % away from midline (outside), to midline (inside)
                    % y: intrusion, -y: extrusion
                    % +z: distal to 1-1, 2-1, -z: mesal to 1-8, 2-8
                    % determine vector to next tooth in order
                    for g=1:2
                        gidx=find(toothGroupV==g);
                        [~,tidx]=sort(toothNumberV(gidx),'descend');
                        sidx=gidx(tidx);
                        toothVectorV(:,sidx(1:end-1)) = toothCenterV(:,sidx(2:end))-toothCenterV(:,sidx(1:end-1));
                    end
                    % do front teeth, from 1-1 to 2-1, from 2-1 to 1-1
                    gidx=find(toothNumberV==1);
                    [~,tidx]=sort(toothGroupV(gidx),'descend');
                    sidx=gidx(tidx);
                    toothVectorV(:,sidx(1:end-1)) = toothCenterV(:,sidx(2:end))-toothCenterV(:,sidx(1:end-1));
                    [~,tidx]=sort(toothGroupV(gidx),'ascend');
                    sidx=gidx(tidx);
                    toothVectorV(:,sidx(1:end-1)) = toothCenterV(:,sidx(2:end))-toothCenterV(:,sidx(1:end-1));
                else
                    % x: away from midline (outside), to midline  lingual (inside)
                    % y: intrusion; -y: extrusion
                    % +z: distal to 3-n -> 3-(n-1), 4-1, -z: mesal to 3-8, 4-8
                    for g=3:4
                        gidx=find(toothGroupV==g);
                        [~,tidx]=sort(toothNumberV(gidx),'descend');
                        sidx=gidx(tidx);
                        toothVectorV(:,sidx(1:end-1)) = toothCenterV(:,sidx(2:end))-toothCenterV(:,sidx(1:end-1));
                    end
                    % do two front teeth
                    sidx=find(toothNumberV==1);
                    toothVectorV(:,sidx(1)) = toothCenterV(:,sidx(2))-toothCenterV(:,sidx(1));
                    toothVectorV(:,sidx(2)) = toothCenterV(:,sidx(1))-toothCenterV(:,sidx(2));
                end
                % show vectors to next tooth
                figure(329)
                hold on
                quiver3(toothCenterV(1,:),toothCenterV(2,:),toothCenterV(3,:),toothVectorV(1,:),toothVectorV(2,:),toothVectorV(3,:),0)
                hold off
                legend(legStr)
                
                %% setup local tooth coordinate system
                % 3 orthogonal vectors, xz defined by toothVectorV
                % normalized to 1
                maxC = 5;
                aW=maxC/4;  % arrow width
                aH=maxC/2;  % arrow head size
                
                % show gum and teeth
                % -1: show all in cornsilkColor for graphical abstract
                %  0: top view
                %  1: front view
                %  2: only transparent teeth
                if saveGAbstractFigs || saveFigs
                    coorLoopV=-1:2;
                else
                    coorLoopV=0:2;
                end
                if (saveFigs || saveMainFigs) && endN==refN
                    for coorLoop=coorLoopV
                        figure(510+coorLoop),clf
                        if coorLoop<2
                            % gum and teeth
                            fh=trisurf(endTR);
                            if coorLoop==1
                                disp(['Figure ' num2str(510+coorLoop) ': gum and teeth, front view'])
                                % was this
                                set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                                % make it brighter like other view
                                set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',0.7*fAlpha,'FaceLighting',fLighting)
                            elseif coorLoop==-1
                                disp(['Figure ' num2str(510+coorLoop) ': gum and teeth in cornsilkColor for graphical abstract'])
                                % all in cornsilkColor for graphical abstract
                                set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                            else
                                disp(['Figure ' num2str(510+coorLoop) ': gum and teeth, top view'])
                                % was this
                                set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',0.7*fAlpha,'FaceLighting',fLighting)
                                % to make it darker like front view, *darker*
                                set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                            end
                            hold on
                            % show teeth end state
                            fh=trisurf(endTeethTR);
                            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',0.7*fAlpha,'FaceLighting',fLighting)
                            % add wisdom teeth
                            if endN==refN
                                for to=1:length(seg{p}.extraToothStrV)
                                    fh=trisurf(seg{p}.extraEndToothTR{to});
                                    set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',0.7*fAlpha,'FaceLighting',fLighting)
                                end
                            end
                        else
                            disp(['Figure ' num2str(510+coorLoop) ': only teeth'])
                            % only transparent teeth
                            fh=trisurf(endTeethTR);
                            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',0.7*fAlpha,'FaceLighting',fLighting)
                            hold on
                        end
                        if coorLoop>-1
                            % add markers of occlusion plane
                            plot3(occPlane.P1(1),occPlane.P1(2),occPlane.P1(3),'k*','MarkerSize',ms-3,'LineWidth',lw)
                            plot3(occPlane.P2(1),occPlane.P2(2),occPlane.P2(3),'k*','MarkerSize',ms-3,'LineWidth',lw)
                            plot3(occPlane.P3(1),occPlane.P3(2),occPlane.P3(3),'k*','MarkerSize',ms-3,'LineWidth',lw)
                        end
                        
                        hold off
                        axis equal
                        if coorLoop<1
                            % from top
                            if doUpper
                                [caz,cel]=view([0 -1 0]);
                            else
                                [caz,cel]=view([0 -1 0]);
                            end
                        elseif coorLoop==1
                            % from front
                            [caz,cel]=view([0 0 1]);
                        else
                            if doUpper
                                [caz,cel]=view([1 -1 0]);
                            else
                                [caz,cel]=view([-1 -1 0]);
                            end
                        end
                        lightangle(gca,caz,cel)
                        material(materialV)
                        axis off
                        
                        % make it same limits
                        xlim(xlimRot*1.2)
                        ylim(ylimRot-[15 0])
                        zlim(zlimRot*1.2)
                    end
                    
                    % add arrows
                    figure(327), clf
                    figure(326), clf
                    figure(325), clf
                end
                %
                flipXM = eye(3); flipXM(1,1)=-1;  % matrix to flip in x direction
                toothCoordV = zeros(3,3,noT);
                for to=1:noT
                    % in-plane vector
                    tmpV = toothVectorV([1 3],to)/norm(toothVectorV([1 3],to));
                    v3 = [tmpV(1); 0; tmpV(2)];
                    % orthogonal vector in occlusion plane
                    v1 = [-v3(3); 0; v3(1)];
                    % y, stays in occlusion plane
                    v2 = cross(v1,v3);
                    
                    th=-90*pi/180;
                    if doUpper && toothGroupV(to)==2
                        % flip one side to orient lingual->buccal in x
                        toothCoordV(:,:,to)=[v1 v2 v3]*flipXM;
                        th=90*pi/180;
                    elseif doUpper==0 && toothGroupV(to)==3
                        toothCoordV(:,:,to)=[v1 v2 v3]*flipXM;
                        % rotation matrix for buccular axis
                        th=90*pi/180;
                    else
                        toothCoordV(:,:,to)=[v1 v2 v3];
                    end
                    % save local transformations
                    seg{p}.toothCoordV = toothCoordV;
                    
                    % transform individual teeth around their center
                    dX=bsxfun(@minus,seg{p}.endToothTR{to}.Points,toothCenterV(:,to)');   % position same cente
                    newT=dX*toothCoordV(:,:,to)+toothCenterV(:,to)';
                    % transform vector to next tooth
                    % then pointing towards mesal
                    newVec=toothVectorV(:,to)'*toothCoordV(:,:,to);
                    
                    if (saveFigs || saveMainFigs) && endN==refN
                        % show original teeth from top
                        figure(325)
                        hold on
                        fh=trisurf(seg{p}.endToothTR{to});
                        set(fh,'FaceColor',seg{p}.color{to},'EdgeColor','none','FaceAlpha',fAlphaLight,'FaceLighting',fLighting)
                        
                        % use arrow3 function with filled arrowhead
                        % arrow in mesial direction
                        mesVec = maxC*toothCoordV(:,3,to)';
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+mesVec,'k2',maxC/2,maxC)
                        % arrow in buccal direction
                        buccVec = maxC*toothCoordV(:,1,to)';
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+buccVec,'k2',maxC/2,maxC)
                        % arrow in extrusion direction, i.e. -intrusion
                        extVec = -maxC*toothCoordV(:,2,to)';
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+extVec,'k2',maxC/2,maxC)
                        hold off
                        
                        for coorLoop=0:2
                            % add local coordinate vectors
                            figure(510+coorLoop)
                            hold on
                            % use arrow3 function with filled arrowhead
                            % arrow in mesial direction, light twany (dark yellow)
                            mesVec = 1.7*maxC*toothCoordV(:,3,to)';
                            arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+mesVec,'^t2',aW,aH)
                            
                            % arrow in buccal direction, blue
                            buccVec = 1.7*maxC*toothCoordV(:,1,to)';
                            arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+buccVec,'b2',aW,aH)
                            
                            % arrow in extrusion direction, dark violet
                            extVec = -1.7*maxC*toothCoordV(:,2,to)';
                            arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+extVec,'_v2',aW,aH)
                            
                            if to==1
                                % add coordinate axis to have same arrows
                                occYVec = 1.7*maxC*[0 0 1];
                                occYVecLen = norm(abs(occPlane.P1-occPlane.midP));
                                arrow3(occPlane.midP,occPlane.midP+[0 0 1]*occYVecLen+occYVec,'k2',aW,aH)
                                % arrow in x direction
                                occXVec = 1.7*maxC*[1 0 0];
                                if doUpper
                                    arrow3(occPlane.P2,occPlane.P3+occXVec*0.8,'k2',aW,aH)
                                else
                                    if doXFlipLowerDisp
                                        % fix arrow
                                        arrow3(occPlane.P3,occPlane.P2+occXVec*0.8,'k2',aW,aH)
                                    else
                                        arrow3(occPlane.P2,occPlane.P3-occXVec*0.8,'k2',aW,aH)
                                    end
                                end
                                % arrow in extrusion direction, i.e. -intrusion
                                xlim(occPlane.midP(1)+xlimMidP)
                                ylim(occPlane.midP(2)+ylimMidP)
                                if p==2
                                    zlim(occPlane.midP(3)+zlimMidP)
                                else
                                    zlim(occPlane.midP(3)+zlimMidP-5)  % to fit widsom tooth
                                end
                            end
                            hold off
                        end
                        
                        figure(325)
                        if to==noT
                            axis equal tight
                            if doUpper
                                [caz,cel]=view([0 -1 0]);
                            else
                                [caz,cel]=view([0 1 0]);
                            end
                            lightangle(gca,caz,cel)
                            if makeClean
                                axis off
                            else
                                xlabel('x')
                                ylabel('y')
                                zlabel('z')
                                grid on
                                title('original teeth')
                                text(toothCenterV(1,to),2*toothCenterV(2,to),toothCenterV(3,to),num2str(seg{p}.toothStrV{to}))
                            end
                        end
                        
                        % show original teeth from side
                        figure(326)
                        hold on
                        fh=trisurf(seg{p}.endToothTR{to});
                        set(fh,'FaceColor',seg{p}.color{to},'EdgeColor','none','FaceAlpha',fAlphaLight,'FaceLighting',fLighting)
                        
                        % use arrow3 function with filled arrowhead
                        % arrow in mesial direction
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+mesVec,'k2',maxC/2,maxC)
                        % arrow in buccal direction
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+buccVec,'k2',maxC/2,maxC)
                        % arrow in extrusion directio
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+extVec,'k2',maxC/2,maxC)
                        hold off
                        if to==noT
                            axis equal tight
                            [caz,cel]=view([1 -1 0]);
                            lightangle(gca,caz,cel)
                            if makeClean
                                axis off
                            else
                                xlabel('x')
                                ylabel('y')
                                zlabel('z')
                                grid on
                                title('original teeth')
                                text(toothCenterV(1,to),2*toothCenterV(2,to),toothCenterV(3,to),num2str(seg{p}.toothStrV{to}))
                            end
                        end
                        
                        % show anatomically oriented teeth
                        figure(327)
                        hold on
                        newTR = triangulation(seg{p}.endToothTR{to}.ConnectivityList,newT);
                        fh=trisurf(newTR);
                        set(fh,'FaceColor',seg{p}.color{to},'EdgeColor','none','FaceAlpha',fAlphaLight-0.1,'FaceLighting',fLighting)
                        % use arrow3 function with filled arrowhead
                        newVec = maxC*newVec/norm(newVec);
                        arrow3(toothCenterV(:,to)',toothCenterV(:,to)'+newVec,'k2',maxC/2,maxC)
                        hold off
                        if to==noT
                            axis equal tight
                            xlabel('buccal->')
                            ylabel('intrusion->')
                            zlabel('mesial->')
                            grid on
                            % view onto occlusion xz-plane
                            if doUpper
                                [caz,cel]=view([0 -1 0]);
                            else
                                [caz,cel]=view([0 1 0]);
                            end
                            lightangle(gca,caz,cel)
                            title('anatomically oriented teeth')
                            text(toothCenterV(1,to),2*toothCenterV(2,to),toothCenterV(3,to),num2str(seg{p}.toothStrV{to}))
                        end
                    end
                end
                
                if saveFigs && endN==refN
                    if makeClean
                        print(325,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_xyzCS_clean.png'])
                        print(326,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_xyzCS_side_clean.png'])
                        print(327,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_anatomyCS_clean.png'])
                    else
                        print(325,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_xyzCS.png'])
                        print(326,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_xyzCS_side.png'])
                        print(327,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg_anatomyCS.png'])
                    end
                end
                
                if saveFigs && endN==refN
                    % for graphical abstract
                    if doUpper==0 && doXFlipLowerDisp
                        % flip in x and rotate 180, as flip in y
                        print(509,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        % flipud
                        imwrite(flipud(toXFlipI),[regFigsDir pat{p}.id '_' partStr '_' endStr '_clean.png'])
                        % shines through if lighter!
                        figure(509)
                        % set background color,
                        set(gcf,'color',gumColor)
                        set(509, 'InvertHardCopy', 'off');
                        print(509,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        % flipud
                        imwrite(flipud(toXFlipI),[regFigsDir pat{p}.id '_' partStr '_' endStr '_gumBG_clean.png'])
                    else
                        print(509,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(toXFlipI,[regFigsDir pat{p}.id '_' partStr '_' endStr '_clean.png'])
                        % shines through if lighter!
                        figure(509)
                        % set background color,
                        set(gcf,'color',gumColor)
                        set(509, 'InvertHardCopy', 'off');
                        print(509,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        % no flipping
                        imwrite(toXFlipI,[regFigsDir pat{p}.id '_' partStr '_' endStr '_gumBG_clean.png'])
                    end
                end
                
                if (saveFigs || saveMainFigs) && endN==refN
                    if doUpper==0 && doXFlipLowerDisp
                        % Fig 1 in manuscript
                        % flip in x and rotate 180, as flip in y
                        print(510,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(flipud(toXFlipI),[regFigsDir 'Fig1_' pat{p}.id '_' partStr '_' endStr '_teeth_gum_xyzCS_darker.png'])
                        
                        print(512,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(flipud(toXFlipI),[regFigsDir 'Fig1_' pat{p}.id '_' partStr '_' endStr '_teeth_xyzCS.png'])
                        % flip in x and rotate 180, as flip in y
                        print(511,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(flipud(toXFlipI),[regFigsDir 'Fig1_' pat{p}.id '_' partStr '_' endStr '_teeth_gum_front_xyzCS_brighter.png'],'BitDepth',8)
                    else
                        % no up-down flipping
                        print(510,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(toXFlipI,[regFigsDir 'Fig1_' pat{p}.id '_' partStr '_' endStr '_teeth_gum_xyzCS_darker.png'])
                        print(512,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(toXFlipI,[regFigsDir 'Fig1_' pat{p}.id '_' partStr '_' endStr '_teeth_xyzCS.png'])
                        print(511,'-dpng','tmpFlipX.png','-r300')
                        toXFlipI=imread('tmpFlipX.png');
                        imwrite(toXFlipI,[regFigsDir 'Fig1_' pat{p}.id '_' partStr '_' endStr '_teeth_gum_front_xyzCS_brighter.png'],'BitDepth',8)
                    end
                end
                %%
            end
            if ~runThrough
                keyboard
            end
            %break
            
            %% visualize individual teeeth
            disp('Figure 330: show individual teeth')
            figure(330)
            fh=trisurf(endTeethTR);
            set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
            legendStr{1}='Teeth';
            hold on
            for i=1:length(seg{p}.endToothTR)
                fh=trisurf(seg{p}.endToothTR{i});
                set(fh,'FaceColor',seg{p}.color{i},'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                legendStr{i+1}=[seg{p}.toothStrV{i} ' ' seg{p}.short{i}];
            end
            hold off
            axis equal
            % stick to fixed image axis
            xlim(xlimSave)
            ylim(ylimSave)
            zlim(zlimSave)
            
            xlabel('x')
            ylabel('y')
            zlabel('z')
            [caz,cel]=view([1 -1 0]);
            lightangle(gca,caz,cel)
            material(materialV)
            title(endStr)
            legend(legendStr,'Location','bestoutside')
            if saveFigs
                print(330,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' endStr '_teethSeg.png'])
            end
            if ~runThrough
                keyboard
            end
            
            %% determine motion and error per tooth
            for to=1:noT
                % find closest points of given tooth in overall mesh as
                % non-rigid results are for all teeth
                [idxV,distV] = knnsearch(endTeethTR.Points, seg{p}.endToothTR{to}.Points);
                if max(distV)>0
                    % should not happen
                    disp(['max distance ' num2str(max(distV))])
                    keyboard
                end
                
                % show in endTR state with rigid(stateTRhat) aligned
                % invNrTransfPointsFull: non-rigid transform endTR to
                % stateTR to get same meshes (stateTRhat)
                invNrTransfPointsTooth = invNrTransfPointsFull(idxV,:);
                nrEnd2StateToothTRglobalCS=triangulation(seg{p}.endToothTR{to}.ConnectivityList,invNrTransfPointsTooth);
                % rigid transformation of stateTRhat to endTR
                rigidStateToothPoints=transformPointsForward(tform,nrEnd2StateToothTRglobalCS.Points);
                seg{p}.rigidState2EndToothTRglobalCS{to} = triangulation(nrEnd2StateToothTRglobalCS.ConnectivityList,rigidStateToothPoints);
                if alignTeethAnatomy
                    % transform endTR to anatomy coordinate system
                    dX=bsxfun(@minus,seg{p}.endToothTR{to}.Points,toothCenterV(:,to)');
                    endTRPointsACS=dX*seg{p}.toothCoordV(:,:,to)+toothCenterV(:,to)';
                    seg{p}.endToothTRanatomyCS{to}=triangulation(seg{p}.endToothTR{to}.ConnectivityList,endTRPointsACS);
                    % transform rigid(stateTRhat) to anatomy coordinate system
                    dX=bsxfun(@minus,rigidStateToothPoints,toothCenterV(:,to)');
                    rigidStateToothPoints=dX*seg{p}.toothCoordV(:,:,to)+toothCenterV(:,to)';
                    seg{p}.rigidState2EndToothTRanatomyCS{to} = triangulation(nrEnd2StateToothTRglobalCS.ConnectivityList,rigidStateToothPoints);
                end
                if doEndSimReg
                    % nrMovDisp: from simTR to rigid(endTR), i.e.
                    % non-rigid(endTR) - rigid(endTR)
                    seg{p}.stateToothTRglobalCS{to}=triangulation(seg{p}.endToothTR{to}.ConnectivityList,...
                        seg{p}.endToothTR{to}.Points-nrMovDisp(idxV,:));
                else
                    % nrMovDisp: from stateTR to rigid(endTR), i.e.
                    % rigid(endTR) - non-rigid(EndTR)
                    seg{p}.stateToothTRglobalCS{to}=triangulation(seg{p}.endToothTR{to}.ConnectivityList,...
                        seg{p}.endToothTR{to}.Points-nrMovDisp(idxV,:));
                end
                
                %% show initial and final state in global coordCS
                % for manuscript, for endN=10
                % for tooth 11, 13, 15, 16, 31, 33, 35, 36
                showCentered = 1;
                if (saveFigs || saveMainFigs) && (~doEndSimReg && endN==refN)
                    for dispTLoop=1:2
                        endTMean=mean(seg{p}.endToothTR{to}.Points);
                        disp('Figure 950: show individual teeth overlayed')
                        figure(950)
                        if dispTLoop==1
                            % in anatomyCS
                            tmpToothTR=seg{p}.endToothTRanatomyCS{to};
                        else
                            tmpToothTR=seg{p}.endToothTR{to};
                        end
                        if showCentered
                            % show centered
                            tmpPoints=tmpToothTR.Points;
                            for i1=1:3
                                tmpPoints(:,i1)=tmpPoints(:,i1)-endTMean(i1);
                            end
                            fh=trisurf(triangulation(tmpToothTR.ConnectivityList,tmpPoints));
                        else
                            fh=trisurf(tmpToothTR);
                        end
                        set(fh,'FaceColor',BMCredLighter,'EdgeColor','none','FaceAlpha',fAlphaLight-0.2,'FaceLighting',fLighting)
                        hold on
                        if dispTLoop==1
                            % in anatomyCS
                            tmpToothTR=seg{p}.rigidState2EndToothTRanatomyCS{to};
                        else
                            tmpToothTR=seg{p}.rigidState2EndToothTRglobalCS{to};
                        end
                        if showCentered
                            % show centered
                            tmpPoints=tmpToothTR.Points;
                            for i1=1:3
                                tmpPoints(:,i1)=tmpPoints(:,i1)-endTMean(i1);
                            end
                            fh=trisurf(triangulation(tmpToothTR.ConnectivityList,tmpPoints));
                        else
                            fh=trisurf(tmpToothTR);
                        end
                        set(fh,'FaceColor',BMCblueLighter,'EdgeColor','none','FaceAlpha',fAlphaLight-0.2,'FaceLighting',fLighting)
                        hold off
                        % show intrusion and buccular
                        if doUpper
                            [caz,cel]=view([0 0 1]);
                        else
                            [caz,cel]=view([0 0 -1]);
                        end
                        axis off equal
                        % same size around mean of endTR
                        if showCentered
                            xlim([-7 7])
                            ylim([-7 7])
                            zlim([-7 7])
                            set(gca,'XTick',[-5 0 5])
                            set(gca,'YTick',[-5 0 5])
                            set(gca,'ZTick',[-5 0 5])
                        else
                            xlim([-7 7]+endTMean(1))
                            ylim([-7 7]+endTMean(2))
                            zlim([-7 7]+endTMean(3))
                        end
                        if saveFigs
                            if dispTLoop==1
                                printFname = [regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_viewZ_clean.png'];
                            else
                                printFname = [regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_viewZ_clean.png'];
                            end
                            print(950,'-dpng',printFname)
                        end
                        pause(1)
                        axis on
                        if dispTLoop==1
                            % in anatomy
                            xlabel('buccal\rightarrow')
                            if doUpper
                                ylabel('\leftarrowextrusion')
                            else
                                ylabel('extrusion\rightarrow')
                            end
                            zlabel('mesial\rightarrow')
                            % color axes, same color as arrows from local coordinate
                            % tooth system
                            ax=gca;
                            ax.XColor='blue';
                            arrow3Violet = [0.6191 0 0.6191];
                            ax.YColor=arrow3Violet;
                            arrow3Tawny = [0.8747 0.5235 0];
                            ax.ZColor=arrow3Tawny;
                            ax.Box='on';
                            ax.BoxStyle='full';
                        else
                            xlabel('x [mm]')
                            ylabel('y [mm]')
                            zlabel('z [mm]')
                        end
                        set(gca,'FontSize',fs)
                        if saveFigs
                            if dispTLoop==1
                                printFname = [regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_viewZ.png'];
                            else
                                printFname = [regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr  '_teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_viewZ.png'];
                            end
                            drawnow
                            print(950,'-dpng',printFname)
                        end
                        if saveFigs || saveMainFigs
                            if dispTLoop==1
                                % show in 3D with y-axis up
                                [caz,cel]=view(3);
                                cel=90-cel;
                                if doUpper
                                    view(caz,cel)
                                    camup([0 1 0])
                                    ylabel('\leftarrowextrusion')
                                    xlabel('buccal\rightarrow','Rotation',-5)
                                    zlabel('mesial\rightarrow','Rotation',-caz+10)
                                else
                                    caz=180+caz;
                                    view(caz,cel)
                                    camup([0 -1 0])
                                    ylabel('extrusion\rightarrow')
                                    zlabel('mesial\rightarrow','Rotation',180-caz+10)
                                    
                                    if buccLowerToRight
                                        set(gca,'xdir','reverse')
                                    else
                                        xlabel('\leftarrowbuccal','Rotation',-5)
                                    end
                                end
                                lightangle(gca,caz,cel)
                                material(materialV)
                                set(gca,'units','normalized','position',[0 0 1 1]);
                                
                                if doUpper
                                    printFname = [regFigsDir 'Fig4_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both.png'];
                                else
                                    if buccLowerToRight
                                        printFname = [regFigsDir 'Fig5_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_buccToRight.png'];
                                    else
                                        printFname = [regFigsDir 'Fig5_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both.png'];
                                    end
                                end
                                drawnow
                                print(950,'-dpng',printFname)
                                axis off
                                if doUpper
                                    printFname = [regFigsDir 'Fig4_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_clean.png'];
                                else
                                    if buccLowerToRight
                                        printFname = [regFigsDir 'Fig5_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_buccToRight_clean.png'];
                                    else
                                        printFname = [regFigsDir 'Fig5_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_clean.png'];
                                    end
                                end
                                drawnow
                                print(950,'-dpng',printFname)
                                
                                set(fh,'BackFaceLighting','unlit')
                                if doUpper
                                    printFname = [regFigsDir 'Fig4_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_unlitback_clean.png'];
                                else
                                    if buccLowerToRight
                                        printFname = [regFigsDir 'Fig5_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_unlitback_buccToRight_clean.png'];
                                    else
                                        printFname = [regFigsDir 'Fig5_' pat{p}.id '_' partStr '_' stateStr 'To' endStr '_' alignStr 'teeth_rigidReg_tooth' seg{p}.toothStrV{to} '_both_unlitback_clean.png'];
                                    end
                                end
                                drawnow
                                print(950,'-dpng',printFname)
                            end
                        end
                    end
                    if ~runThrough
                        keyboard
                    end
                end
                
                %% fit a local rigid transformation
                % to motion of tooth from non-rigid registration of all teeth
                % by
                % - calculating centroids
                % - aligning centroids with translation
                % - then rotating
                % transforms A to fit B, i.e. transfA2B = R*A+t
                % center of rotation should be tooth center in endTR!!!
                
                % in endTR state with rigid(stateTR) to fit endTR
                if alignTeethAnatomy
                    B=seg{p}.endToothTRanatomyCS{to}.Points'-toothCenterV(:,to);
                    A=seg{p}.rigidState2EndToothTRanatomyCS{to}.Points'-toothCenterV(:,to);
                else
                    B=seg{p}.endToothTR{to}.Points'-toothCenterV(:,to);
                    A=seg{p}.rigidState2EndToothTRglobalCS{to}.Points'-toothCenterV(:,to);
                end
                if doEndSimReg
                    % want endTR -> simN transformation, hence swap
                    Acopy=A;
                    A=B;
                    B=Acopy;
                end
                
                % determine 3D rigid transformation
                [seg{p}.Rot{to}, seg{p}.transl{to}, lrms] = Kabsch(A,B);
                
                tmpM=logm(seg{p}.Rot{to});
                % rotations around x, y, z
                % http://mesh.brown.edu/rotations/
                rotDegV=360*[tmpM(3,2) tmpM(1,3) tmpM(2,1)]/(2*pi);
                
                % check application of found rigid transformation
                % transformation acts like this
                transfA2B=seg{p}.Rot{to}*A+seg{p}.transl{to};
                meanErr=mean(sum((B-transfA2B).^2,1));
                sqrtErr=sqrt(meanErr);  % should be the same as lrms!
                if abs(sqrtErr-lrms)>eps
                    disp(['A2B error ' num2str(sqrtErr,'%.6f') ' vs. ' num2str(lrms,'%.6f') ' mm'])
                end
                transfA2BTR = triangulation(seg{p}.endToothTR{to}.ConnectivityList,transfA2B'); % local rigid transf fitting stateTR
                mainRotDeg = norm(rotDegV);  mainTransl = norm(seg{p}.transl{to});
                % output found local rigid transformation, LRMS fitting error, residual
                if to==1
                    disp('& ID    & [deg] & [deg] & [deg] & [deg] & [mm] & [mm] & [mm] & [mm] & [mm]')
                    if alignTeethAnatomy
                        disp('& tooth &  buc  &  int  &  mes  & main  & buc  & int  & mes   & main & RMS')
                    else
                        disp('& tooth &  rx   &  ry   &  rz   & main  & tx   & ty   & tz   & main & RMS')
                    end
                end
                disp([tableStr ' & ' seg{p}.toothStrV{to} ' ' num2str(rotDegV,' & %.2f') ' ' num2str(mainRotDeg,' & %.2f') ' ' num2str(seg{p}.transl{to}',' & %.2f') ' ' num2str(mainTransl,' & %.2f') ' & ' num2str(lrms,'%.3f') ' \\'])
                
                %% show displacements between center of teeth
                if doEndSimReg==1 && length(doEndSimRegV)>1
                    % calculate (relative) achieved tooth motion
                    % displacement vectors at tooth center
                    
                    % directly as in figure 810+to in endTR state
                    if alignTeethAnatomy
                        % in anatomical coordinate system
                        P1=mean(saveStartTR.seg.rigidState2EndToothTRanatomyCS{to}.Points); %State_n
                        P2=mean(seg{p}.endToothTRanatomyCS{to}.Points);  % endTR
                        P3=mean(seg{p}.rigidState2EndToothTRanatomyCS{to}.Points);  % Sim_n
                    else
                        % in occlusal plane xyz coordinate system
                        P1=mean(saveStartTR.seg.rigidState2EndToothTRglobalCS{to}.Points); %State_n
                        P2=mean(seg{p}.endToothTR{to}.Points);  % endTR
                        P3=mean(seg{p}.rigidState2EndToothTRglobalCS{to}.Points);  % Sim_n
                    end
                    dx1=P2-P1;   % state_n->End
                    dx2=P3-P2;   % End->SimEnd
                    dx3=dx1+dx2; % state_n->SimEnd
                    
                    % calculate relative achievement
                    normV(1)=norm(dx1);
                    normV(2)=norm(dx2);
                    normV(3)=norm(dx3);
                    % remaining to do, End->Sim_n relative to State_n->Sim_n
                    achievedPerc=100*(1-normV(2)/normV(3));
                        
                    % determine transformation parameters
                    % centered rigid transformation, with center at endTR (toothCenterV(:,to))
                    % achieved T1 (State_n->End)
                    T1=[saveStartTR.seg.Rot{to} saveStartTR.seg.transl{to}; 0 0 0 1]; 
                    tmpMT1=logm(saveStartTR.seg.Rot{to});
                    startRotDegV=180*[tmpMT1(3,2) tmpMT1(1,3) tmpMT1(2,1)]/pi;
                    % missing T2(End->Sim_endN), 
                    T2=[seg{p}.Rot{to} seg{p}.transl{to}; 0 0 0 1]; 
                    tmpM=logm(seg{p}.Rot{to});
                    rotDegV=180*[tmpM(3,2) tmpM(1,3) tmpM(2,1)]/pi;
                    % planned T3 (State_n->End->Sim_endN)
                    T3=T2*T1;
                    tmpMT3=logm(T3(1:3,1:3));
                    rotDegT3=180*[tmpMT3(3,2) tmpMT3(1,3) tmpMT3(2,1)]/pi;
                    translT3=T3(1:3,4)';
                    
                    % difference missing: planned - achieved
                    
                    % achieved T1 (State_n->End)
                    seg{p}.achievedDegM(:,to)=startRotDegV;
                    seg{p}.achievedNormDegV(to)=norm(startRotDegV);
                    seg{p}.achievedTransM(:,to)=saveStartTR.seg.transl{to}';
                    seg{p}.achievedNormTransV(to)=norm(saveStartTR.seg.transl{to});
                    % missing T2(End->Sim_endN)
                    seg{p}.missingDegM(:,to)=rotDegV;
                    seg{p}.missingNormDegV(to)=norm(rotDegV);
                    seg{p}.missingTransM(:,to)=seg{p}.transl{to}';
                    seg{p}.missingNormTransV(to)=norm(seg{p}.transl{to});
                    
                    % planned motion T3 (State_n->Sim_endN)
                    seg{p}.plannedDegM(:,to)=rotDegT3;
                    seg{p}.plannedNormDegV(to)=norm(rotDegT3);
                    seg{p}.plannedTransM(:,to)=translT3;
                    seg{p}.plannedNormTransV(to)=norm(translT3);
                       
                    mainRotDeg = norm(rotDegV);  mainTransl = norm(seg{p}.transl{to});
                       
                    disp([pat{p}.id ' T1 Start' num2str(n) '$\rightarrow$' endStr ' & ' seg{p}.toothStrV{to} ' ' num2str(startRotDegV,' & %.2f') ' ' num2str(saveStartTR.seg.transl{to}',' & %.2f') ' \\'])
                    disp([pat{p}.id ' T2 ' tableStr ' & ' seg{p}.toothStrV{to} ' ' num2str(rotDegV,' & %.2f') ' ' num2str(seg{p}.transl{to}',' & %.2f') ' \\'])
                    disp([pat{p}.id ' T3 Start' num2str(n) '$\rightarrow$Sim' num2str(endN) '  & ' seg{p}.toothStrV{to} ' ' num2str(rotDegT3,' & %.2f') ' ' num2str(translT3,' & %.2f') ' \\'])
                    disp([pat{p}.id ' Achieved  & ' seg{p}.toothStrV{to} ' ' num2str(seg{p}.achievedDegM(:,to)',' & %.2f') ' & ' num2str(seg{p}.achievedNormDegV(to),'%.2f') ' ' num2str(seg{p}.achievedTransM(:,to)',' & %.2f') ' & ' num2str(seg{p}.achievedNormTransV(to),'%.2f') ' \\'])
                    disp([pat{p}.id ' Missing   & ' seg{p}.toothStrV{to} ' ' num2str(seg{p}.missingDegM(:,to)',' & %.2f') ' & ' num2str(seg{p}.missingNormDegV(to),'%.2f') ' ' num2str(seg{p}.missingTransM(:,to)',' & %.2f') ' & ' num2str(seg{p}.missingNormTransV(to),'%.2f') ' \\'])
                    disp([pat{p}.id ' Planned   & ' seg{p}.toothStrV{to} ' ' num2str(seg{p}.plannedDegM(:,to)',' & %.2f') ' & ' num2str(seg{p}.plannedNormDegV(to),'%.2f') ' ' num2str(seg{p}.plannedTransM(:,to)',' & %.2f') ' & ' num2str(seg{p}.plannedNormTransV(to),'%.2f') ' \\'])
                end
            end
            if ~runThrough
                keyboard
            end
            %% save results of stateTR to calculate relative error
            if doEndSimReg==0
                saveStartTR.seg = seg{p};
            end
            
            if n~=refN
                % show propagated teeth on state
                figure(51),clf
                fh=trisurf(rigidMovingTR);
                if doUpper==0
                    rotate(fh,[1 0 0],180);
                end
                set(fh,'FaceColor',gumColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                hold on
                % show teeth in state space
                for to=1:noT
                    fh=trisurf(seg{p}.stateToothTRglobalCS{to});
                    if doUpper==0
                        rotate(fh,[1 0 0],180);
                    end
                    set(fh,'FaceColor',cornsilkColor,'EdgeColor','none','FaceAlpha',fAlpha,'FaceLighting',fLighting)
                end
                hold off
                axis equal
                % from top
                if doUpper
                    [caz,cel]=view([0 -1 0]);
                else
                    [caz,cel]=view([0 1 0]);
                end
                lightangle(gca,caz,cel)
                set(gcf,'color',lightBMCblue)
                set(51, 'InvertHardCopy', 'off');
                material(materialV)
                axis off
                
                % make it same limits
                xlim(xlimRot)
                ylim(ylimRot)
                zlim(zlimRot)
                
                if saveFigs
                    print(51,'-dpng',[regFigsDir pat{p}.id '_' partStr '_' stateStr 'To' endStr '_rigidReg_gum.png'])
                end
            end
            
        end
    end   % loop stateTR, endSim
    if ~runThrough
        keyboard
    end
    %% combined statistics over same categories
    if isfield(seg{p},'plannedDegM')
        
        noCat=length(toothCat);
        % statistics of planned motion
        resPlanned = NaN(noCat,8);
        resNumbers = zeros(noCat,1);
        for cat=1:noCat
            fidx=find(strcmp(seg{p}.name,toothCat{cat}.name));
            meanDegCat = mean(abs(seg{p}.plannedDegM(:,fidx)),2,'omitnan');
            meanNormDegCat = mean(abs(seg{p}.plannedNormDegV(fidx)),2,'omitnan');
            meanTranslCat = mean(abs(seg{p}.plannedTransM(:,fidx)),2,'omitnan');
            meanNormTransCat = mean(abs(seg{p}.plannedNormTransV(fidx)),2,'omitnan');
            resPlanned(cat,:) = [meanDegCat' meanNormDegCat' meanTranslCat' meanNormTransCat'];
            resNumbers(cat) = length(fidx);
            disp([num2str(resNumbers(cat)) ' & '  toothCat{cat}.name ' ' num2str(meanDegCat',' & %.2f') ...
                ' & ' num2str(meanNormDegCat,' %.2f') ...
                ' ' num2str(meanTranslCat',' & %.2f') ' & ' num2str(meanNormTransCat,' %.2f') ' \\'  ])
        end
        % total
        disp([num2str(sum(resNumbers)) ' & meanPlanned ' num2str(mean(resPlanned,'omitnan'),' & %.2f') ' \\'  ])
        
        
        % statistics of missing motion, i.e. absolute difference
        disp('absolute missing')
        resMissing = NaN(noCat,8);
        resNumbers = zeros(noCat,1);
        for cat=1:noCat
            fidx=find(strcmp(seg{p}.name,toothCat{cat}.name));
            meanDegCat = mean(abs(seg{p}.missingDegM(:,fidx)),2,'omitnan');
            meanNormDegCat = mean(abs(seg{p}.missingNormDegV(fidx)),2,'omitnan');
            meanTranslCat = mean(abs(seg{p}.missingTransM(:,fidx)),2,'omitnan');
            meanNormTransCat = mean(abs(seg{p}.missingNormTransV(fidx)),2,'omitnan');
            resMissing(cat,:) = [meanDegCat' meanNormDegCat' meanTranslCat' meanNormTransCat'];
            resNumbers(cat) = length(fidx);
            disp([num2str(resNumbers(cat)) ' & '  toothCat{cat}.name ' ' num2str(meanDegCat',' & %.2f') ...
                ' & ' num2str(meanNormDegCat,' %.2f')...
                ' ' num2str(meanTranslCat',' & %.2f') ' & ' num2str(meanNormTransCat,' %.2f') ' \\'  ])
        end
        % total
        disp([num2str(sum(resNumbers)) ' & MeanMissing ' num2str(mean(resMissing,'omitnan'),' & %.2f') ' \\'  ])
        
        if saveData
            % save teeth errors, achieved, planned, missing
            saveFname = [resultDir pat{p}.id '_' partStr '_State' num2str(n) stateStr 'To' endStr '_absTeethResults.mat'];
            patseg = seg{p};
            disp(['save results to ' saveFname])
            save(saveFname,'patseg')
        end
    end
    %%
    if ~runThrough
        keyboard
    end
end % loop lower, upper jaw
end % for state n=1-9 and endN=N; or n=1 and endN=2:N
end % loop patients p

