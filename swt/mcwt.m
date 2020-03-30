% MCWT �����߶Ȳ�������morletС���ֽ�����Щ�߶��µ�����С��ϵ��
%   coefs  - С���任ϵ������ÿ�ж�Ӧһ���߶ȵķֽ�ϵ������Ӧ��Ƶ���ν���
%   data   - �����ź�,���ݵ�����Ϊż��
%   wtype  - ����С�������ͣ�'real'����'complex'
%   scales - �źŵķֽ�߶ȣ���Ϊ��һ�߶Ȼ��߶�
%   dt     - ��������
%   ds     - С�������ĳ߶Ȳ������
%
% author - Wei Wang
% �ο����ף�Christopher Torrence and Gilbert P.Compo. A practical Guide to wavelet Analysis[J]. Bulletin of American Meteorological Society, 1998,79(1): 61-78.

function coefs = mcwt(data,scales,wave)

N = length(data);  % �źŲ�����������Ϊż��  data�ǲ���������ݣ�N�ǲ���������ݳ��ȣ���MCA-bcr��
J = length(scales); % С���ֽ�ĳ߶��� ����������ȷ����

% �źŵ�Ƶ���ʾ
Fdata = zeros(N,1);
Fdata = fft(data(:)); % �ź������ٸ���Ҷ�任 fft

coefs = zeros(N,J); % С���任��ϵ������ ��ʼ��Ϊ�����
am=1.0;
% С���任������
for j = 1:J
    % ����߶�Ϊscales(j)��С������ wave()��ÿ���߶��µ�morletС���ĸ���Ҷ�任 ��ô�õ��ģ������������� 
    coefs(:,j) = ifft(Fdata.*conj(wave(:,j))).*am;   %% conj()������ĺ���   �������  Ƶ������� Ȼ������任
end


