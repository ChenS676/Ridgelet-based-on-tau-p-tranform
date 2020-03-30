%  WAVELET 
% this function give the discrete morlet function 
% Reference: % �ο����ף�Christopher Torrence and Gilbert P.Compo. A practical Guide to
% wavelet Analysis[J]. Bulletin of American Meteorological Society, 1998,79(1): 61-78.
%  3. Wavelet analysis b. Wavelet transform paragraph 1



function wave=waveLet(N,dt,scales,J)
% ���������źŵ�Ƶ����
% dw = 2*pi/dt/N;     % Ƶ��������
% fpl = fix(N/2);     % Ƶ�������᳤��   ��0��£ȡ��
% fnl = fix((N-1)/2); % Ƶ�򸺰��᳤��
% k  = (1:fpl).*dw;
% 
% w1 = [0, k];       % ������
% w1 = w1(:);
% w2 = -k(fnl:-1:1); % ������
% w2 = w2(:);
% am = 0.5;  %%%%%%%%%amʲô��˼?����ʲô�ã�
% 
% % morletС������
% w0 = 6;            % ����Ƶ��   where ��0 is the nondimensional frequency,here taken to be 6 to satisfy the admissibility condition


w0 = 6;%����Ƶ��
dw = 2*pi/dt/N;%Ƶ��������
w1 = [0:N/2]*dw;
w2 = [N/2+1:N-1]*dw;
wave  = zeros(N,J); % Ƶ��morletС������

% С���任������
for j = 1:J
    wave(1:N/2+1,j) = sqrt(2*pi*scales(j)/dt)*pi^(-0.25)*exp(-(scales(j)*w1-w0).^2./2);
    wave(N/2+2:N,j) = sqrt(2*pi*scales(j)/dt)*pi^(-0.25)*exp(-(-scales(j)*w2-w0).^2./2);
end