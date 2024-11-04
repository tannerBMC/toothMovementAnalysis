% data definition and registration parameters for 
% the aligner performance case study for two patients reported in 
% 
% I. Filippon, C. Tanner, J. A. von Jackowski, G. Schulz, T. Toepper, B.
% Mueller, University of Basel, 2024
% Aligner-induced tooth movements in three dimensions using clinical data of two patients
%
% imaging, tooth and  patient data file definitions
% definition for 2 patients 
%
% May 2024, Christine Tanner, Biomaterials Science Center, University of Basel

%% figure definitions

BMCred =[198 19  24]/255;    % hex #C61318, complement is green [14 152 46]; 3 colors: #D16BAS (pink), #5FFBF1 (
BMCblue=[0   77 145]/255;    % hex #014D91, complement is orange [178 78 6
BMCmag=(BMCred+BMCblue)/2;   %
								      
lightBMCblue=BMCblue/4+[0.75 0.75 0.75];  % add white
            
BMCgreen=[14 152 46]/255;
BMCorange=[178 78 6]/255;
BMCyellow=(BMCgreen+BMCorange)/2;   %
% make lighter
hsv1=rgb2hsv(BMCred);
hsv1(3)=hsv1(3)*1.1; % 10% lighter
BMCredLighter=hsv2rgb(hsv1);
%
hsv1=rgb2hsv(BMCblue);
hsv1(3)=hsv1(3)*1.1; %
BMCblueLighter=hsv2rgb(hsv1);
%
cornsilkColor  = [255 248 220]/255;
wheatColor     = [255 222 173]/255;
gumColor = [0.5 0.5 0.5]+([255 87 51]/510);  % white and pink
             
% color map with 8 colors
cmap=prism(8);
colorV = {'y','r','m','b','c','k','g','y','r','m','b','c','k','g'};
colorLineV = {'y-','r-','m-','b-','c-','k-','g-','y--','r--','m--','b--','c--','k--','g--'};
colorMarkersV = {'yo','rx','m*','bd','cv','k^','g+','yo','rx','m*','bd','cv','k^','g+'};

coordStrV = {'x','y','z'};

% plotting parameters
fs=15;   % FontSize
ms=15;   % MarkerSize
lw=3;    % LineWidth

% face properties
fAlpha = 1.0;
fAlphaLight=0.8;
% calculates the vertex normals and interpolates linearly across the faces. Select this method to view curved surfaces
fLighting = 'gouraud';
% for less saturated specles
materialV=[0.4 0.4 0.4];

%% define data

baseDir = ['data' filesep];
figsDir = [baseDir 'figs' filesep];

patV = {'3485','6457'};
noP=length(patV);

refN=10;   % final state T9 of first nine weeks with weekly intraoral scans, used as reference
N=11;      % total number of teeth states, T0:T9 weekly and T10 for final state after many weeks
stateStr = cell(1,N);
for n=1:N
    stateStr{n}=['T' num2str(n-1)];
end

% for subjects 1 & 2
% CutM: define manually which isoline should be used for cutting teeth
% 0 unknown yet, stops to show user possiblities
% -1 dont cut (contours not closed!), read instead manually cut mesh _cut1
% StraightV: 0: PCA defined cut angle, 1: keep straight
% KeepHighV: 0: keep below threshold, 1: keep above threshold

%
% LOWER - mandible
%
% default is cut1, i.e. set -1 for cut1
stateLowerCutM = -1*ones(N,noP);
stateSimLowerCutM = -1*ones(N,noP);

% default, 0: PCA defined cut angle
stateLowerStraightM = zeros(N,noP);
stateSimLowerStraightM = zeros(N,noP);

% default, 0: keep below threshold, 1: keep above threshold
stateLowerKeepHighM = zeros(N,noP);
stateSimLowerKeepHighM = ones(N,noP);
					      
% largest body mesh
stateLowerCleanM = zeros(N,noP);
stateSimLowerCleanM = zeros(N,noP);

%
% UPPER - maxillary
%
% default, set -1 for cut1
stateUpperCutM = -1*ones(N,noP);
stateSimUpperCutM = -1*ones(N,noP);

% default, 0: PCA defined cut angle
stateUpperStraightM = zeros(N,noP);
stateSimUpperStraightM = zeros(N,noP);
% set 1 for keep straight
stateSimUpperStraightM([1 5],1) = 1;
stateSimUpperStraightM(:,2) = 1;

% default, 0: keep below threshold, 1: keep above threshold
stateUpperKeepHighM = zeros(N,noP);
stateSimUpperKeepHighM = ones(N,noP);
					      
% largest body mesh
stateUpperCleanM = zeros(N,noP);
stateSimUpperCleanM = zeros(N,noP);

% default no tooth segmentations are available-1
% initialize patient definitions
pat = cell(noP,1);
seg = repmat({struct('upperToothStrV',[],'lowerToothStrV',[],'toothStrV',[])}, noP, 1);

%% define occlusion planes by 3 points for endTR
% P1: between front teeth
% P2: distobuccular hoecker on 2nd molar on anatomical right side (upper 1-7, lower 4-7) 
% P3: distobuccular hoecker on 2nd molar on anatomical left side  (upper 2-7, lower 3-7) 

% clear all
upperOccPlane = cell(1,noP);
lowerOccPlane = cell(1,noP);

% set manually determined points [x y z]
upperOccPlane{1}.P1 = [-9.8421   -8.2143   26.0685]; 
upperOccPlane{1}.P2 = [-27.8510   -6.6035  -15.4140]; 
upperOccPlane{1}.P3 = [25.8475   -3.8557   -2.5125]; 

lowerOccPlane{1}.P1 = [-7.3285   -7.6276   25.0078]; 
lowerOccPlane{1}.P2 = [26.3526 -9.18039 -1.72466];  
lowerOccPlane{1}.P3 = [-24.716 -6.17749 -13.3343];  

upperOccPlane{2}.P1 = [-1.9318   -9.1141   26.3856]; 
upperOccPlane{2}.P2 = [-32.5699  -10.3990  -13.6267]; 
upperOccPlane{2}.P3 = [30.3429  -10.2828  -14.3173]; 

lowerOccPlane{2}.P1 = [-1.4049  -10.1022   20.1666]; 
lowerOccPlane{2}.P2 = [29.6122   -6.1890  -17.0080]; 
lowerOccPlane{2}.P3 = [-28.1837   -8.2212  -19.2044]; 

%% define available tooth segmentations per patient

% teeth segmentations for last time step
seg{1}.upperToothStrV = {'17','16','15','13','12','11','21','22','23','25','26','27'};
seg{1}.lowerToothStrV = {'37','36','35','33','32','31','41','42','43','45','46','47'};

% wisdom teeth only to show 
seg{1}.extraUpperToothStrV = {'18','28'};
seg{1}.extraLowerToothStrV = {'48'};

seg{2}.upperToothStrV = {'17','16','15','14','13','12','11','21','22','23','24','25','26','27'};
seg{2}.lowerToothStrV = {'37','36','35','34','33','32','31','41','42','43','44','45','46','47'};
seg{2}.extraUpperToothStrV = {};
seg{2}.extraLowerToothStrV = {};

% tooth categories and colors
colorHSV=hsv(8);
toothCat{1}.name='Central Incisor';   
toothCat{1}.short='CI';
toothCat{1}.color = colorHSV(1,:); 
toothCat{2}.name='Lateral Incisor';   
toothCat{2}.short='LI';
toothCat{2}.color = colorHSV(2,:); 
toothCat{3}.name='Canine';            
toothCat{3}.short='Ca';
toothCat{3}.color = colorHSV(3,:); 
toothCat{4}.name='Premolar1';         
toothCat{4}.short='PM1';
toothCat{4}.color = colorHSV(4,:); 
toothCat{5}.name='Premolar2';        
toothCat{5}.short='PM2';
toothCat{5}.color = colorHSV(5,:); 
toothCat{6}.name='Molar1';           
toothCat{6}.short='M1';
toothCat{6}.color = colorHSV(6,:); 
toothCat{7}.name='Molar2';           
toothCat{7}.short='M2';
toothCat{7}.color = colorHSV(7,:);
   
% do not include wisdom teeth in evaluation             
toothCat{8}.name='Molar3';            
toothCat{8}.short='M3';
toothCat{8}.color = colorHSV(8,:); 

%% define scan names per patient

for p=1:noP
    pat{p}.id = patV{p};
    pat{p}.dir = [pat{p}.id '_Clinical_Trial' filesep 'Sirona' filesep];
    pat{p}.simDir = [pat{p}.id '_Clinical_Trial' filesep 'Bottmedical' filesep];
    % for registration
    %
    % upper - maxillary
    %
    for n=1:N
        % define all 10 states
        if stateUpperCutM(n,p)==-1
            pat{p}.stateUpper{n} = [stateStr{n} filesep pat{p}.id '_OnyxCeph3_Export_OK_cut1'];
            pat{p}.stateUpperOrig{n} = [stateStr{n} filesep pat{p}.id '_OnyxCeph3_Export_OK'];
        else
            pat{p}.stateUpper{n} = [stateStr{n} filesep pat{p}.id '_OnyxCeph3_Export_OK'];
        end
        % clean up, one body only
        if stateUpperCleanM(n,p)
            pat{p}.stateUpper{n} = [pat{p}.stateUpper{n} '_1B'];
        end
        % names of propagated teeth segmentation
        pat{p}.stateUpperTeeth{n} = [stateStr{n} '/' pat{p}.id '_OK'];
        pat{p}.stateUpperToothStr{n} = [stateStr{n} '/' pat{p}.id '_Z'];
        pat{p}.stateUpperTeethFine{n} = [stateStr{n} '/' pat{p}.id 'fine_OK'];
        pat{p}.stateUpperToothFineStr{n} = [stateStr{n} '/' pat{p}.id 'fine_Z'];
    end
    % teeth segmentation only for end state
    pat{p}.endUpperTeeth = [stateStr{refN} '/teethSeg/' pat{p}.id '_OK'];
    pat{p}.endUpperToothStr = [stateStr{refN} '/teethSeg/' pat{p}.id '_Z'];
    %
    % lower - mandibular - unterkiefer (UK)
    %
    for n=1:N
        % define all N states
        if stateLowerCutM(n,p)==-1
            pat{p}.stateLower{n} = [stateStr{n} filesep pat{p}.id '_OnyxCeph3_Export_UK_cut1'];
            pat{p}.stateLowerOrig{n} = [stateStr{n} filesep pat{p}.id '_OnyxCeph3_Export_UK'];
        else
            pat{p}.stateLower{n} = [stateStr{n} filesep pat{p}.id '_OnyxCeph3_Export_UK'];
        end
        if stateLowerCleanM(n,p)
            pat{p}.stateLower{n} = [pat{p}.stateLower{n} '_1B'];
        end
        % names of propagated teeth segmentation
        pat{p}.stateLowerTeeth{n} = [stateStr{n} '/' pat{p}.id '_UK'];
        pat{p}.stateLowerToothStr{n} = [stateStr{n} '/' pat{p}.id '_Z'];
    end
    % teeth only defined for end state
    pat{p}.endLowerTeeth = [stateStr{refN} filesep 'teethSeg' filesep pat{p}.id '_UK'];
    pat{p}.endLowerToothStr = [stateStr{refN} filesep 'teethSeg' filesep pat{p}.id '_Z'];
    %
    % simulated PLAN
    %
    for n=1:N
        % define all N states
        if stateSimUpperCutM(n,p)==-1
            pat{p}.stateSimUpper{n} = [stateStr{n} filesep 'Upper_Aligner_' stateStr{n}(2:end) '_cut1'];
            pat{p}.stateSimUpperOrig{n} = [stateStr{n} filesep 'Upper_Aligner_' stateStr{n}(2:end)];
        else
            pat{p}.stateSimUpper{n} = [stateStr{n} filesep 'Upper_Aligner_' stateStr{n}(2:end)];
        end
        if stateSimUpperCleanM(n,p)
            pat{p}.stateSimUpper{n} = [pat{p}.stateSimUpper{n} '_1B'];
        end
        if stateSimLowerCutM(n,p)==-1
            pat{p}.stateSimLower{n} = [stateStr{n} filesep 'Lower_Aligner_' stateStr{n}(2:end) '_cut1'];
            pat{p}.stateSimLowerOrig{n} = [stateStr{n} filesep 'Lower_Aligner_' stateStr{n}(2:end)];
        else
            pat{p}.stateSimLower{n} = [stateStr{n} filesep 'Lower_Aligner_' stateStr{n}(2:end)];
        end
        if stateSimLowerCleanM(n,p)
            pat{p}.stateSimLower{n} = [pat{p}.stateSimLower{n} '_1B'];
        end
    end
    				      
    % define needed rotations for
    % initialize to zero rotation, N states x measured/simulated x upper/lower 
    % rotation defined in degrees, converted to rad later
    pat{p}.rotXM = zeros(N,2,2);
    pat{p}.rotYM = zeros(N,2,2);
    pat{p}.rotZM = zeros(N,2,2);
    if p==1
        % patient 1
        % upper
        pat{p}.rotXM(1:N-3,1,1) = 90;  % measured   % default
        pat{p}.rotXM(N-2,1,1) = 180;  

        pat{p}.rotXM(N-1,1,1) = 90;

        pat{p}.rotXM(1:N,2,1) = -90;          % simulated default
        
        pat{p}.rotXM(11,1,1) = 0;
        pat{p}.rotYM(11,1,1) =60;
        pat{p}.rotZM(11,1,1) = -80;
        
        % lower
        pat{p}.rotXM(1:N-3,1,2) = 90;
        pat{p}.rotXM(N-1:N,1,2) = 90;

        pat{p}.rotYM(11,1,2) = -90;
       
        pat{p}.rotXM(1:N,2,2) = 90;   % simulated default
    elseif p==2
        % patient 2
        % upper
        pat{p}.rotXM(1:N-3,1,1) = 90; % measured default
        pat{p}.rotXM(5,1,1) = -100;
        pat{p}.rotYM(5,1,1) = 180;
        pat{p}.rotZM(5,1,1) = 40;
        %
        pat{p}.rotYM(9,1,1) = 90;
        pat{p}.rotZM(9,1,1) = 180;
        pat{p}.rotXM(N-1:N,1,1) = 90;
        
        pat{p}.rotXM(1:N,2,1) = -90;          % simulated default
        
        % lower
        pat{p}.rotXM(1:N-3,1,2) = 90;  % default
        pat{p}.rotXM(5,1,2) = 45;
        pat{p}.rotYM(5,1,2) = -90;
        pat{p}.rotXM(N-1:N,1,2) = 90;
       
        pat{p}.rotXM(1:N,2,2) = 90;   % simulated default
    else
        disp(['unknown patient ' num2str(p)])
        break
    end
    
end
%% define limits of figures from previous experience

% this covers all
xlimSave=[-37 32];
ylimSave=[-13 6];
zlimSave=[-30 30];

% similar range as with full
xlimTeethSave=[-35 35];
ylimTeethSave=[-10 5];
zlimTeethSave=[-15 30];  % cut as fewer teeth

% limits around occPlane.midP
xlimMidP = [-39 39];
ylimMidP = [-10 15];
zlimMidP = [-10 50];

%% define mesh alignment
% 
% make fixed mesh upright such that it can be easily cut
% rotate around X2 to minimize extend of points

% align end state with occlusion plane
alignOccPlane = 1;     % 1: align end state with occlusion plane
alignTeethAnatomy = 1; % 1: anatomically orient each tooth in end state
if ~alignOccPlane
    alignTeethAnatomy = 0;  % cannot do this without aligning occlusion plane
end
if alignTeethAnatomy
    alignStr = 'anatomyCS_';
else
    alignStr = '';
end
if alignOccPlane
    % y-axis aligned with normal of occlusion plane
    % x-axis aligned with P2-P3
    configStr = 'occPX_';
    extraOutStr = 'cut1, occPX';
else
    configStr = '';
    extraOutStr = 'cut1';
end

%% encode registration parameters

% rigid ICP parameters
doSmallerRigidTol = 1;
defaultTol = [0.01 0.5];
if doSmallerRigidTol
    tolVal = defaultTol/16;  % convergence tolerance
else
    tolVal = defaultTol/4;  % convergence tolerance
end
inlierRatio = 0.9;      % ratio of inlier
maxIter = 500;          % maximum number of iterations
downPercentage = 0.5;

%% registration output file name

regParStr = ['cut1_' configStr 'down' num2str(downPercentage) '_inlier' num2str(inlierRatio) '_maxIter' num2str(maxIter) '_tol' num2str(tolVal(1)) '_' num2str(tolVal(2))];
regParStr(regParStr=='.') = 'p';
    
% non-rigid registration options

%   Options : structured object with fields:
%       gamm : real valued, weights differences in the rotational and skew
%           part of the deformation against the translational part.
%           default 1
%       epsilon : real values, tolerence for change in transformation.
%           default 1e-4
%       lambda : If using the bi-directional distance metric this weights
%           the contribution of the target -> source term.
%           default 1
%       alphaSet : decreasing vector of real-valued stiffness parameters.
%           High stiffness parameters force global transformations whereas
%           low values allow for local deformations.
%           default linspace(100,10,20)
%       biDirectional : logical, specifies that a bi-directional distance
%           is used.
%           default 0
%       useNormals : logical, specifies that surface normals are to be used
%           to project the source onto the target surface. If this term is
%           used then the Source input should contain a normals field.
%           default 0
%       noFinalSnapping : if useNormals==0, then default is that points are
%           snapped to the closest point on the target surface. Setting
%           this option keeps the points regularized as during the
%           registration, i.e. no snapping.
%           default 0
%       plot : logical, specifies that the transformations should be
%           plotted.
%           default 0
%       rigidInit : logical, specifies that rigid ICP should be performed
%           first before allowing non-rigid and non-global deformations.
%           default 1
%       ignoreBoundary, default 1
%       normalWeighting, default 1

% Specify that surface normals are available and can be used.
% If useNormals == 1, then projects transformed points after registration
% along surface normals to target surface, otherwise snap to closest
% points on target. Takes long.
%
nrOptions.useNormals = 0;
% Specify that the source deformations should be plotted.
nrOptions.plot = 1;
% zero weight if normals have not similar angles
if nrOptions.useNormals
    nrOptions.normalWeighting = 1;
else
    nrOptions.normalWeighting = 0;
    nrOptions.noFinalSnapping = 1; % keep regularized
end
if nrOptions.noFinalSnapping
    % allow more deformations
    nrOptions.epsilon = 1;
    nrOptions.alphaSet = linspace(100,10,6);
    nrRegParStr = 'noSnap';   % only 3 refinements
else
    % converge faster
    nrOptions.epsilon = 1;
    nrOptions.alphaSet = linspace(100,50,3);
    nrRegParStr = 'quick';   % only 3 refinements
end

nrOptions.rigidInit = 0;   % use own ICP result

if nrOptions.useNormals
    nrRegParStr = [nrRegParStr '_withNormals'];
end

regParStr = [regParStr '_' nrRegParStr];


