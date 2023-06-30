% Author: Xiang Shen(shen@apm.ac.cn)
function [cPath, cFold, cFile, cName, nFile] = getFilesInFolder(mainDir, sPre, sPost)

if isempty(sPre)
    sPre = '*';
end
if strcmpi(sPost,'*')
    sPost = '';
end

nFile = 0;
BUFFER = 1E6;
cPath = cell(BUFFER,1);
cFold = cell(BUFFER,1);
cFile = cell(BUFFER,1);
cName = cell(BUFFER,1);

cPost = strsplit(sPost,'|'); 
for p=1:length(cPost)
    sPostPart = cPost{p};
    cSubFolds = getSubfoldersInFolder(mainDir, {}); 
    for i=1:length(cSubFolds)
        subdirpath = [cSubFolds{i}, ['\*',sPostPart] ];
        structFile = dir(subdirpath); 
        for j=1:length(structFile)
            sFold = structFile(j).folder;
            sFile = structFile(j).name;
            
            if (structFile(j).isdir==0) && (sPre(1) == '*' || strncmpi(sPre, sFile, length(sPre)))
                nFile = nFile+1;
                cPath{nFile} = [sFold, '\', sFile];
                cFold{nFile} = [sFold, '\'];
                cFile{nFile} = sFile;
                
                sName = sFile(1:end-length(sPostPart));
                if strcmpi(sName(end),'.') 
                    sName = sName(1:end-1);
                end
                cName{nFile} = sName;
            end
        end
    end
end

b = cellfun('isempty',cPath);
cPath = cPath(~b);
cFold = cFold(~b);
cFile = cFile(~b);
cName = cName(~b);
nFile = length(cName);

end