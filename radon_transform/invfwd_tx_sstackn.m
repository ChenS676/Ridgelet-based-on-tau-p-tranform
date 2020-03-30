function b=invfwd_tx_sstackn(m,dt,p,x,type)  

if strcmp(type,'linear')
    calculate='tz(j)+p(i)*x(k)';               %叠加路径为直线
elseif strcmp(type,'hyper')
    calculate='sqrt(tz(j)^2+p(i)^2*x(k)^2)';   %叠加路径为双曲线
elseif strcmp(type,'parab')
    calculate='tz(j)+p(i).^2*x(k).^2';         %叠加路径为抛物线
end

nt=size(m,1);
np=length(p);
nx=length(x);
tz=0*dt:dt:(nt-1)*dt;                       %时间点
data=zeros(nt,nx);
for i=1:np
    for j=1:nt
        for k=1:nx        
          t=eval(calculate);
          tc=t/dt;
          itc=floor(tc);
          it=itc+1;
          fra=tc-itc;
          if it>0&&it<nt
                 data(it,k)=data(it,k)+(1.-fra)*m(j,i);
                 data(it+1,k)=data(it+1,k)+fra*m(j,i);
          end
        end
    end
end
b=data;             