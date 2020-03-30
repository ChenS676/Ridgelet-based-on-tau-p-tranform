%  WAVELET 
% this function give the discrete morlet function 
% Reference: % 参考文献：Christopher Torrence and Gilbert P.Compo. A practical Guide to
% wavelet Analysis[J]. Bulletin of American Meteorological Society, 1998,79(1): 61-78.
%  3. Wavelet analysis b. Wavelet transform paragraph 1



function wave=waveLet(N,dt,scales,J)
% 建立分析信号的频率轴
% dw = 2*pi/dt/N;     % 频域采样间隔
% fpl = fix(N/2);     % 频域正半轴长度   向0靠拢取整
% fnl = fix((N-1)/2); % 频域负半轴长度
% k  = (1:fpl).*dw;
% 
% w1 = [0, k];       % 正半轴
% w1 = w1(:);
% w2 = -k(fnl:-1:1); % 负半轴
% w2 = w2(:);
% am = 0.5;  %%%%%%%%%am什么意思?会有什么用？
% 
% % morlet小波参数
% w0 = 6;            % 中心频率   where ω0 is the nondimensional frequency,here taken to be 6 to satisfy the admissibility condition


w0 = 6;%中心频率
dw = 2*pi/dt/N;%频域采样间隔
w1 = [0:N/2]*dw;
w2 = [N/2+1:N-1]*dw;
wave  = zeros(N,J); % 频域morlet小波函数

% 小波变换主程序
for j = 1:J
    wave(1:N/2+1,j) = sqrt(2*pi*scales(j)/dt)*pi^(-0.25)*exp(-(scales(j)*w1-w0).^2./2);
    wave(N/2+2:N,j) = sqrt(2*pi*scales(j)/dt)*pi^(-0.25)*exp(-(-scales(j)*w2-w0).^2./2);
end