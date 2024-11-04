function colorM = ownColors(colM,colP,doLogSpace)
% create linearly interpolated colormaps
% based on definining base colors and their position
%
% colorM = ownColors(colM,startList,doLogSpace)
%
% colM       color matrix, [Cx3], defining starting colors
% colP       color position [Cx1]. in [0,255], defining color position
% doLogSpace 1: logarithmic change of color [default linear]
%
% colorM   colormap definition [256x3]
%
N=256;
F=1.2;
colorM=zeros(N,3);
if nargin<3
    doLogSpace=0;
end

for c=1:length(colP)-1
    startC = colM(c,:);
    endC = colM(c+1,:);
    % linearly interpolate
    lenC = colP(c+1)-colP(c)+1;
    for i=1:3
        if doLogSpace
            % factor F makes it steeper
            expV=linspace(exp(startC(i))*F,exp(endC(i))*F,lenC);
            colorM(colP(c):colP(c+1),i)=log(expV)/F;
        else
            colorM(colP(c):colP(c+1),i)=linspace(startC(i),endC(i),lenC);
        end
    end
end
if colP(c+1)<N
    % extend last color if not at last position
    lenC = N-colP(c+1)+1;
    for i=1:3
        colorM(colP(c+1):N,i)=endC(i);
    end
end