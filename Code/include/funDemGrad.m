% Author: Xiang Shen(shen@apm.ac.cn)
function [rGradX, rGradY] =  funDemGrad(rDem, Ref)
nX = size(rDem,2);
nY = size(rDem,1);
rGradX = (rDem(:,3:nX)-rDem(:,1:(nX-2)))/(2*Ref.CellExtentInWorldY);
rGradY = (rDem(1:(nY-2),:)-rDem(3:nY,:))/(2*Ref.CellExtentInWorldY);
rGradX = [rGradX(:,1), rGradX, rGradX(:,end)];
rGradY = [rGradY(1,:); rGradY; rGradY(end,:)];
end