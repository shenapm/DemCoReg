% Author: Xiang Shen(shen@apm.ac.cn)
function funPlotDemDif(rDH,sFoldOt,sDemPair,iAlgo,col)

figure;
set(gcf,'Units','inches');
set(gcf,'Position',[2 2 3.5 2.625]);

h = imagesc(rDH);
set(h, 'AlphaData', ~isnan(rDH));
colormap(col);

h = colorbar; 
set(get(h,'label'),'string','Elevation difference/m','FontSize',9);
set(gca, 'CLim', [-30 30]);

axis image; axis off;
sNameOt = [sDemPair,'_', num2str(iAlgo)];
sPathOt = [sFoldOt,sNameOt,'.tif'];
print(sPathOt,'-dtiff','-r600')
end