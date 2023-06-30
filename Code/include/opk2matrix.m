% Author: Xiang Shen(shen@apm.ac.cn)
function R = opk2matrix(opk)
o = opk(1);
p = opk(2);
k = opk(3);

so = sin(o); co = cos(o);
sp = sin(p); cp = cos(p);
sk = sin(k); ck = cos(k);

Ro = [1 0 0; 0 co -so; 0 so co];
Rp = [cp 0 sp; 0 1 0; -sp 0 cp];
Rk = [ck -sk 0; sk ck 0; 0 0 1];

R = Ro * Rp * Rk;
end

