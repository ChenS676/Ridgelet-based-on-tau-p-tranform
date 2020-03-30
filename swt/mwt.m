function [wcoef,scales,nscales,ds] = mwt(s,dt,nt);   %%%小波变换
s0 = 2*dt;                                 % the smallest resolvable scale
ds = 0.125;                                 % scale index sample period, <0.6
J0 = 1;
n=100;
J = log(n*dt/s0)/log(2)/ds;                    % decomposition scales
J = floor(J);
scales = zeros(1,J);
nscales = length(scales);

for jj = J0:J
    scales(jj) = s0*2^((jj-1)*ds);
end
wave = zeros(nt,J);
wave = waveLet(nt,dt,scales,J);
wcoef = mcwt(s,scales,wave);   %%%小波变换
end
