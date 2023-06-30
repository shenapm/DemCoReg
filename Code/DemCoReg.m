% Author:   Tao Li; Xiang Shen(shen@apm.ac.cn)
% Citation: https://tc.copernicus.org/preprints/tc-2022-205/
ProjDef_DemCoReg;

[cPathDemM,~,~,cNameDemM] = getFilesInFolder([sFoldData,'Master','/'],'','.tif');
sNameDemM = cNameDemM{1};
sPathDemM = cPathDemM{1};
sPathAngA = replace(sPathDemM,'.tif','_Aspect.tif');
sPathAngS = replace(sPathDemM,'.tif','_Slope.tif');

[cPathDemS,~,~,cNameDemS] = getFilesInFolder([sFoldData,'Slave','/'],'','.tif');
for i=1:length(cPathDemS)
    sPathDemS = cPathDemS{i};
    sNameDemS = cNameDemS{i};
    
    sDemPair = [sNameDemM,'  ',sNameDemS];
    sPathLog = [sFoldRes, sDemPair, '.log'];
    diary(sPathLog);
    diary on;
    disp('---------------------------------------------------------------------');
    disp(datetime(now,'ConvertFrom','datenum'));
    fprintf('Project: %s\n',sPrj);
    fprintf('Master  Slave: %s\n',sDemPair);
    
    [rDemM, RefDemM] = funImgRead(sPathDemM,-9000,9000);
    [rDemS, RefDemS] = funImgRead(sPathDemS,-9000,9000);
    rAngS = funImgRead(sPathAngS,0,89);
    rAngA = funImgRead(sPathAngA,-9000,9000);
    rDemM(isnan(rAngS+rAngA)) = nan;
    
    %% Pre-processing
    [ixM,iyM] = meshgrid(1:size(rDemM,2), 1:size(rDemM,1));
    [rXM,rYM] = intrinsicToWorld(RefDemM,ixM,iyM);
    rZM = rDemM;
    
    [ixS,iyS] = meshgrid(1:size(rDemS,2), 1:size(rDemS,1));
    [rXS,rYS] = intrinsicToWorld(RefDemS,ixS,iyS);
    rZS = rDemS;
    
    % Image-centered coordinates
    mXS = mean(rXS(:),'omitnan'); rXSc = rXS-mXS; rXMc = rXM-mXS;
    mYS = mean(rYS(:),'omitnan'); rYSc = rYS-mYS; rYMc = rYM-mYS;
    mZS = mean(rZS(:),'omitnan'); rZSc = rZS-mZS; rZMc = rZM-mZS;
    
    % Partitioning
    rng('default');
    nPixel = numel(rDemM);
    cvo = cvpartition(nPixel,'Holdout',0.9);
    bTr = training(cvo); %Training set
    bTe = test(cvo);     %Test set
    
    %% Co-Registration
    for iAlgo = 1:length(cMethod)
        fprintf('Algo #%d %s\n',iAlgo, cMethod{iAlgo});
        
        vMAD = [];
        vabc = [];  %Translation vectors (Cylindrical coordinates)
        vxyz = [];  %Translation vectors (Cartesian coordinates)
        v7pm = [];  %Seven parameters (Scale, rotation, and, translation)
        
        for iter=1:tIterMaxNum
            fprintf('\tIter #%d  ',iter);
            
            if iter==1
                rDemI = rDemS; % Intermediate data
                RefDemI = RefDemS;
                rDemW = MapWarp(rDemI,RefDemI,RefDemM,'linear');
                
                [ixI,iyI] = meshgrid(1:size(rDemI,2), 1:size(rDemI,1));
                [rXI,rYI] = intrinsicToWorld(RefDemI,ixI,iyI);
                rZI = rDemI;
                
                % Image-centered coordinates
                rXIc = rXI-mXS;
                rYIc = rYI-mYS;
                rZIc = rZI-mZS;
                
                % Scattered points
                pIc = [rXIc(:), rYIc(:), rZIc(:)];
                pIc = pIc(~isnan(pIc(:,3)),:);
            else % Updating
                if iAlgo<=3
                    RefDemI.XWorldLimits = RefDemI.XWorldLimits + vxyz(end,1);
                    RefDemI.YWorldLimits = RefDemI.YWorldLimits + vxyz(end,2);
                    rDemI = rDemI+vxyz(end,3);
                    rDemW = MapWarp(rDemI,RefDemI,RefDemM,'linear');
                else
                    R = opk2matrix(v7pm(end,2:4)); % Rotation matrix
                    pIc = (1+v7pm(end,1))*pIc*R' + repmat(v7pm(end,5:7),size(pIc,1),1);
                    pI = pIc + repmat([mXS, mYS, mZS],size(pIc,1),1);
                    
                    [ixW,iyW] = worldToIntrinsic(RefDemM,pI(:,1),pI(:,2));
                    pZW = griddata(ixW,iyW,pI(:,3),ixM,iyM);
                    rDemW = reshape(pZW,size(rDemM,1),size(rDemM,2));
                    rDemW(isnan(rAngS+rAngA+rDemM))=nan;
                end
            end
            
           %% Height differences
            rDH = rDemM-rDemW;
           
            dh = rDH(~isnan(rDH+rAngS+rAngA));
            mdh = median(dh);
            sdh = 1.4826*mad(dh,1);
            dh = dh(abs(dh-mdh)<=3*sdh);
            rDH(abs(rDH-mdh)>3*sdh)=nan;
            vMAD = [vMAD; median(abs(dh))];
            fprintf('MedAD = %7.3f  ',vMAD(end));
            
            if iter==1
                if iAlgo == 1
                    funPlotDemDif(-rDH,sFoldRes,sDemPair,0,colDiv);
                    close all;
                end
            else
                if abs(vMAD(end)-vMAD(end-1)) < tIterMadDif
                    funPlotDemDif(-rDH,sFoldRes,sDemPair,iAlgo,colDiv);
                    close all;
                    break;
                end
            end
            
           %% Preparing training data
            bV = ~isnan(rDH+rAngS+rAngA) & rAngS>=tMinAngS; %Valid pixels
            bTrV = bTr & bV(:); %Valid training samples
            
            dh = rDH(bTrV);
            aspect = rAngA(bTrV);
            slope  = rAngS(bTrV);
            ts = tand(slope);
            dhts = dh./ts;
            
            mdh = median(dh);
            mts = median(tand(slope));
            
            % Image-centered coordinates
            Xc = rXMc(bTrV);
            Yc = rYMc(bTrV);
            Zc = rZMc(bTrV);
            
            % Regression coefficients
            vX = sind(aspect).*ts;
            vY = cosd(aspect).*ts;
            vG = vX.*Xc + vY.*Yc + Zc;  %Gamma
            vO = Yc-vY.*Zc;             %Omega
            vP = vX.*Zc-Xc;             %Phi
            vK = vY.*Xc-vX.*Yc;         %Kappa
            
            %% Regresion
            if iAlgo == 1
                fr = nlinfit([aspect,ts],dh,Nuth2,[0 0 mdh],opts);
                a = fr(1); b = fr(2); c = fr(3);
                vabc = [vabc; a, b, c];
                vxyz = [vxyz; a*sind(b), a*cosd(b), c];                
            elseif iAlgo == 2
                fr = nlinfit(aspect,dhts,Nuth3,[0 0 mdh/mts],opts);
                a = fr(1); b = fr(2); c = fr(3)*mts;
                vabc = [vabc; a, b, c];
                vxyz = [vxyz; a*sind(b), a*cosd(b), c];
            elseif iAlgo == 3
                fr = robustfit([vX,vY],dh);
                dX = fr(2); dY = fr(3); dZ = fr(1);
                vxyz = [vxyz; dX, dY, dZ];
                vabc = [vabc; sqrt(dX.^2+dY.^2), atan2d(dX,dY), dZ];
            else
                fr = robustfit([vX,vY,vG,vO,vP,vK],dh);
                dX = fr(2); dY = fr(3); dZ = fr(1);
                vxyz = [vxyz; dX, dY, dZ];
                vabc = [vabc; sqrt(dX.^2+dY.^2), atan2d(dX,dY), dZ];
                v7pm = [v7pm; fr(4),fr(5),fr(6),fr(7),fr(2),fr(3),fr(1)];
            end
            
            % Transformation parameters
            if iAlgo<=3
                fprintf('abc = (%8.3f %8.3f %8.3f) ',vabc(end,1),vabc(end,2),vabc(end,3));
                fprintf('xyz = (%8.3f %8.3f %8.3f) ',vxyz(end,1),vxyz(end,2),vxyz(end,3));
                fprintf('d = %6.3f\n',rms(vxyz(end,:)));
            else
                fprintf('7pm = (%8.1e %8.1e %8.1e %8.1e %8.3f % 8.3f %8.3f)\n',...
                    v7pm(end,1),v7pm(end,2),v7pm(end,3),v7pm(end,4),v7pm(end,5),v7pm(end,6),v7pm(end,7));
            end
        end
        
        % Cumulative values
        if iAlgo<=3
            vxyzA = cumsum(vxyz);
            vxyzA = vxyzA(end,:);
            fprintf('xyzA= (%8.3f %8.3f %8.3f)\n',vxyzA(end,1),vxyzA(end,2),vxyzA(end,3));
        else
            v7pmA = cumsum(v7pm);
            v7pmA = v7pmA(end,:);
            fprintf('7pmA= (%8.1e %8.1e %8.1e %8.1e %8.3f % 8.3f %8.3f)\n',...
                v7pmA(end,1),v7pmA(end,2),v7pmA(end,3),v7pmA(end,4),v7pmA(end,5),v7pmA(end,6),v7pmA(end,7));
        end
        close all;
    end
    diary off;
end