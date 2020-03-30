% MICWT   MORLET INVERSE CONTINOUS WAVELET TRANSFORM 
%   �����߶Ȳ�����С���ֽ�ϵ������morletС���ع���Ӧ���źŷ���
%       Ҫ��߶Ȳ�����˳����С��ϵ�������������һ��
%   yn     - �����߶ȷ�Χ���ع��ź�
%   coefs  - С���任ϵ������ÿ�ж�Ӧһ���߶ȵķֽ�ϵ������Ӧ��Ƶ���ν���
%   scales - �źŵķֽ�߶ȣ���Ϊ��һ�߶Ȼ��߶ȣ���ֽ�����Ӧ
%   dt     - ��������
%   ds     - С�������ĳ߶Ȳ������
% author - Wei Wang
% �ο����ף�Christopher Torrence and Gilbert P.Compo. A practical Guide to
% wavelet Analysis[J]. Bulletin of American Meteorological Society, 1998,
% 79(1): 61-78.

function yn = micwt(coefs,scales,dt,ds)

% �źŲ�������N���ֽ�߶���cJ
% coefs=real(coefs);
[N cJ]= size(coefs);
coefs=real(coefs);
% С���ֽ�ĳ߶���
J = length(scales);

% �ж��������
if cJ~=J
    error('Inconsistent input parameter!');
end

% morletС�������Ĳ���
Cdelta = 0.776;      % addimissibility condition constant
phi0   = pi^(-0.25); % remove energy scaling

% �ź��ع���С����ص�ϵ��
Crec = ds*(dt^0.5)/Cdelta/phi0*0.9968;
% keep in mind that  Cdelta,phi0,Crec  depends the type of the wavelet 

% ��ʼ��ԭʼ�ź�
yn = zeros(N,1);
% С�����任������
for j = 1:J
    % ����߶��ع�ԭʼ�ź�
    % �㷨û����  2017.05.08
    yn = yn + Crec.*real(coefs(:,j)./sqrt(scales(j)));
end

%% 


