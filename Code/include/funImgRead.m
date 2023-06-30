function [r,Ref] = funImgRead(sPath,minv,maxv)
[r, Ref] = readgeoraster(sPath,'OutputType','double');

r(r<minv) = NaN;
r(r>maxv) = NaN;

end
