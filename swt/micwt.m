% MICWT   MORLET INVERSE CONTINOUS WAVELET TRANSFORM 
%   给定尺度参数及小波分解系数，用morlet小波重构对应的信号分量
%       要求尺度参数的顺序与小波系数矩阵的列向量一致
%   yn     - 给定尺度范围的重构信号
%   coefs  - 小波变换系数矩阵，每列对应一个尺度的分解系数，相应主频依次降低
%   scales - 信号的分解尺度，可为单一尺度或多尺度，与分解矩阵对应
%   dt     - 采样周期
%   ds     - 小波函数的尺度采样间隔
% author - Wei Wang
% 参考文献：Christopher Torrence and Gilbert P.Compo. A practical Guide to
% wavelet Analysis[J]. Bulletin of American Meteorological Society, 1998,
% 79(1): 61-78.

function yn = micwt(coefs,scales,dt,ds)

% 信号采样点数N，分解尺度数cJ
% coefs=real(coefs);
[N cJ]= size(coefs);
coefs=real(coefs);
% 小波分解的尺度数
J = length(scales);

% 判断输入参数
if cJ~=J
    error('Inconsistent input parameter!');
end

% morlet小波函数的参数
Cdelta = 0.776;      % addimissibility condition constant
phi0   = pi^(-0.25); % remove energy scaling

% 信号重构与小波相关的系数
Crec = ds*(dt^0.5)/Cdelta/phi0*0.9968;
% keep in mind that  Cdelta,phi0,Crec  depends the type of the wavelet 

% 初始化原始信号
yn = zeros(N,1);
% 小波反变换主程序
for j = 1:J
    % 逐个尺度重构原始信号
    % 算法没看懂  2017.05.08
    yn = yn + Crec.*real(coefs(:,j)./sqrt(scales(j)));
end

%% 


