function [cFold] = getSubfoldersInFolder(mainDir, cFold)

cFold = [cFold; mainDir];

subDir = dir( mainDir );
for i = 1:length( subDir )
    if( isequal( subDir( i ).name, '.' )||...
        isequal( subDir( i ).name, '..')||...
        ~subDir( i ).isdir)  
        continue;
    end
    subFold = fullfile(mainDir, subDir(i).name);
    cFold = getSubfoldersInFolder(subFold, cFold);
end

end