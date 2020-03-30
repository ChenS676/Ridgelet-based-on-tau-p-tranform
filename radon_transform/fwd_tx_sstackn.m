function m=fwd_tx_sstackn(d,dt,p,x,type)

if strcmp(type,'linear')
    calculate='tz(j)+p(i)*x(k)';               %����·��Ϊֱ��
elseif strcmp(type,'hyper')
    calculate='sqrt(tz(j)^2+p(i)^2*x(k)^2)';   %����·��Ϊ˫����
elseif strcmp(type,'parab')
    calculate='tz(j)+p(i).^2*x(k).^2';         %����·��Ϊ������
end

nt=size(d,1);
np=length(p);
nx=length(x);
tz=0:dt:(nt-1)*dt;                         %ʱ���
data=zeros(nt,np);
for i=1:np
    for j=1:nt
        for k=1:nx
            t=eval(calculate);    
            tc=t/dt;
            itc=floor(tc);
            it=itc+1;
            fra=tc-itc;
            if it>0&&it<nt
                 data(j,i)=data(j,i)+(1.-fra)*d(it,k)+fra*d(it+1,k);      %���Բ�ֵ
            end
        end
    end
end
m=data;                 