% Author: Xiang Shen(shen@apm.ac.cn)
function [rS,rA] = funDemSlopeAspect(r,R,sMethod)
nX = size(r,2);
nY = size(r,1);

fX_Z = padarray((r(:,3:nX)-r(:,1:(nX-2)))/(2*R.CellExtentInWorldY),[0,1],nan);
fY_Z = padarray((r(1:(nY-2),:)-r(3:nY,:))/(2*R.CellExtentInWorldY),[1,0],nan);
if strncmpi(sMethod,'Z',1) % ZevenbergenThorne1987
    rS = atand(sqrt(fX_Z.^2+fY_Z.^2));
    rA = 180-atand(fY_Z./fX_Z)+90*fX_Z./abs(fX_Z);
    rA(fX_Z==0) = 180*(fY_Z(fX_Z==0)>0);
    rA(rA>=360) = rA(rA>=360)-360;
else % Horn1981
    fX_H = (fX_Z + movsum(fX_Z,3)  )/4;
    fY_H = (fY_Z + movsum(fY_Z,3,2))/4;
    
    rS = atand(sqrt(fX_H.^2+fY_H.^2));
    rA = 180-atand(fY_H./fX_H)+90*fX_H./abs(fX_H);
    rA(fX_H==0) = 180*(fY_H(fX_H==0)>0);
    rA(rA>=360) = rA(rA>=360)-360;
end
end
