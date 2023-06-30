% Author: Xiang Shen(shen@apm.ac.cn)
function [rW,RW] = MapWarp(rA,RA,RB,method,varargin)

if strcmpi(class(RA),'map.rasterref.MapPostingsReference')
    RA = mapref_postings2cells(RA);
end

if isobject(RB)
    if strcmpi(class(RB),'map.rasterref.MapPostingsReference')
        RB = mapref_postings2cells(RB);
    end
    bbx = [RB.XWorldLimits,RB.YWorldLimits];
else
    bbx = RB;
end

if ~strcmpi(method,'nearest') && ~strcmpi(method,'linear') && ~strcmpi(method,'cubic') && ~strcmpi(method,'spline')
    method = 'linear';
end

if nargin>=5 && ~isempty(varargin{1})
    szPixel = varargin{1};
else
    if isobject(RB)
       szPixel = [RB.CellExtentInWorldX,RB.CellExtentInWorldY];
    else
       szPixel = [RA.CellExtentInWorldX,RA.CellExtentInWorldY];
    end
end

if nargin>=6 && varargin{2}
    minX = floor((bbx(1)-RA.XWorldLimits(1))/szPixel(1))*szPixel(1) + RA.XWorldLimits(1);
    maxX =  ceil((bbx(2)-RA.XWorldLimits(2))/szPixel(1))*szPixel(1) + RA.XWorldLimits(2);
    minY = floor((bbx(3)-RA.YWorldLimits(1))/szPixel(2))*szPixel(2) + RA.YWorldLimits(1);
    maxY =  ceil((bbx(4)-RA.YWorldLimits(2))/szPixel(2))*szPixel(2) + RA.YWorldLimits(2);
    bbx = [minX,maxX,minY,maxY];
end

bbx(2) = bbx(1) + round((bbx(2)-bbx(1))/szPixel(1))*szPixel(1);
bbx(4) = bbx(3) + round((bbx(4)-bbx(3))/szPixel(2))*szPixel(2);

%%
RW = maprefcells(bbx(1:2),bbx(3:4),szPixel(1),szPixel(2),'ColumnsStartFrom','north');
[x,y] = meshgrid(1:RW.RasterSize(2), 1:RW.RasterSize(1));
[X,Y] = intrinsicToWorld(RW,x,y);
b = X>RA.XWorldLimits(1) & X<RA.XWorldLimits(2) & Y>RA.YWorldLimits(1) & Y<RA.YWorldLimits(2);

if size(rA,3)==1
    rW = NaN(RW.RasterSize);
    rW(b) = mapinterp(rA,RA,X(b),Y(b),method);
else
    nBand = size(rA,3);
    rW = NaN([RW.RasterSize,nBand]);
    for bi = 1:nBand
        r = NaN(RW.RasterSize);
        r(b) = mapinterp(rA(:,:,bi),RA,X(b),Y(b),method);
		rW(:,:,bi) = r;
    end
end

funName = str2func(class(rA));
if islogical(rA)
    rW(isnan(rW)) = 0;
end
rW = funName(rW);

end