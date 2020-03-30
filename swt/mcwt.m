% MCWT 给定尺度参数，用morlet小波分解求这些尺度下的连续小波系数
%   coefs  - 小波变换系数矩阵，每列对应一个尺度的分解系数，相应主频依次降低
%   data   - 输入信号,数据点数须为偶数
%   wtype  - 分析小波的类型，'real'或者'complex'
%   scales - 信号的分解尺度，可为单一尺度或多尺度
%   dt     - 采样周期
%   ds     - 小波函数的尺度采样间隔
%
% author - Wei Wang
% 参考文献：Christopher Torrence and Gilbert P.Compo. A practical Guide to wavelet Analysis[J]. Bulletin of American Meteorological Society, 1998,79(1): 61-78.

function coefs = mcwt(data,scales,wave)

N = length(data);  % 信号采样点数，需为偶数  data是补了零的数据，N是补了零的数据长度，在MCA-bcr中
J = length(scales); % 小波分解的尺度数 在主程序中确定了

% 信号的频域表示
Fdata = zeros(N,1);
Fdata = fft(data(:)); % 信号做快速傅里叶变换 fft

coefs = zeros(N,J); % 小波变换的系数矩阵 初始化为零矩阵
am=1.0;
% 小波变换主程序
for j = 1:J
    % 计算尺度为scales(j)的小波函数 wave()是每个尺度下的morlet小波的傅里叶变换 怎么得到的？？？？？？？ 
    coefs(:,j) = ifft(Fdata.*conj(wave(:,j))).*am;   %% conj()是求共轭的函数   卷积定理  频率里相乘 然后做逆变换
end


