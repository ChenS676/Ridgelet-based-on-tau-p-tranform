%% first actual data test
% Processing is as follows
% 1 single channel filter to remove most surface waves % low pass 15Hz
% 2 radon transform cut surface wave residue
% 3 radon The ridgelet transform is used to remove part of the body wave energy and put it back
% 4 last low-pass filtering
addpath 'C:\Users\strohwitwer\Documents\Mein Folder\Ridgelet-based-on-tau-p-transform\swt'
addpath 'C:\Users\strohwitwer\Documents\Mein Folder\Ridgelet-based-on-tau-p-transform\radon_transform'
%% INPUT
close all;clear all; clc;
dt = 0.004; 
nt = 512;
nx = 109;  
dx = 0.005;
t  = (0:nt-1).*dt;
t=t'; 
x = (0:nx-1).*1;

shotground = zeros(nt,nx);
fidin = fopen('datasets\shotground.segy','r');
volume_head = fread(fidin,[3600,1],'*uchar');      % 3600: number of juantou
for i=1:nx
    trace_head = fread(fidin, [240,1], '*uchar');     % 240 number of  shot
    shotground(:,i)=fread(fidin,[nt,1],'float'); 
end
fclose(fidin);
set(0,'defaultfigurecolor','w')
figure,pcolor(1:nx, t, shotground);shading interp; axis ij;colormap('gray'),colorbar
%%  step1:  Preprocess
s = shotground;

close all
for i=1:nx
   fre1 = 15/(1/dt/2);   %低通15Hz
   [B,A] = butter(8,fre1,'low');
   lfils(1:nt,i) = filtfilt(B,A,s(1:nt,i));
end
figure,
pcolor(lfils);shading interp; axis ij;colormap('gray'),colorbar

% output
fr_err = s-lfils;
fbd = fr_err;    

figure,
pcolor(fr_err);shading interp; axis ij;colormap('gray'),colorbar
%% Step 2: RADON TRANSFORM

pmax = 20;  
P=100;
p = 0:pmax/P:pmax;

np = length(p);
Para_1 = struct('dt',dt,'nt',nt,'t',t,'dx',dx,'nx',nx,'x',x,'p',p,'np',np,'pmax',pmax);
type = 'linear'
fra = cgnr_radon(fbd,Para,type);

figure,
pcolor(1:np,1:nt,fra),shading interp;axis ij;colormap('Parula'),colorbar,set(gca,'XAxisLocation','top')
set(gcf, 'Renderer', 'ZBuffer');
%%

bfra = zeros(size(fra));
bfra(:,20:70)= fra(:,20:70);
figure,pcolor(1:np,1:nt,bfra),shading interp;axis ij;colormap('Parula'),colorbar,set(gca,'XAxisLocation','top')
set(gcf, 'Renderer', 'ZBuffer');

%% 重建回时空域
bfragr = invfwd_tx_sstackn(bfra,dt,p,x,type);
figure,imagesc(bfragr),shading interp;axis ij;colormap('gray'),colorbar,set(gca,'XAxisLocation','top'),set(gcf, 'Renderer', 'ZBuffer');

raffbd = fbd - bfragr;
figure,imagesc(raffbd),shading interp;axis ij;colormap('gray'),colorbar,set(gca,'XAxisLocation','top'),set(gcf, 'Renderer', 'ZBuffer');
 subplot(1,3,2),

%% %%   阈值分割 按斜率切片 依照频带分割
close all
for i = 1:P
[mwt_coef,scales,ds] = mwt(bfra(:,i),dt,nt);  % 每个斜率的ridgelet分解系数 mwt_coef 
part_l(:,i) = reshape(mwt_coef,nt*45,1);       % ridgelet分解系数降维
%[gr_reccoef,bd_reccoef,r,c] = thres_l(mwt_coef,160);  
[gr_reccoef,r,c] = thres_l(mwt_coef,'linear',scales); % 线性切割，留下大于40Hz的部分
part_thresgr(:,i) = reshape(gr_reccoef,r*c,1);   % 阈值分割后的信号降维重排
%part_thresbd(:,i) = reshape(bd_reccoef,r*c,1);
end

%% 重建至radon域
rec = zeros(nt,nx);
for j = 1:P
ri_grcoef = reshape(part_thresgr(:,j),r,c);   
ra_grcoef = micwt(ri_grcoef,scales,dt,ds);  % 对可能过度滤波的体波信号进行恢复
ra_gr(:,j) = ra_grcoef;
end 

% input reconstruction 
figure,pcolor(ra_gr);shading interp; axis ij;colormap('gray'),colorbar

% % 重建至信号时空域
%% input rec ,inrec 时空域重建信号
close all
bdrec= invfwd_tx_sstackn(ra_gr,dt,p,x,type);
figure,pcolor(bdrec);shading interp; axis ij;colormap('gray'),colorbar
figure,pcolor(400*bfragr-bdrec);shading interp; axis ij;colormap('gray'),colorbar,%
%% +friraf_bd
figure,pcolor(2000*lfils + bfragr -bdrec);shading interp; axis ij;colormap('gray')
xlabel('Time/s'),ylabel('Trace'),title('Ridgelet低通分离出的体波能量')
set (gca,'position',[0.15,0.15,0.75,0.75],'XAxisLocation','bottom');
set(gcf,'Unit','Normalized','Position',[0.1,0.2,0.38,0.58])
%%
riraf_bd = raffbd + bdrec;
figure,pcolor(riraf_bd);shading interp; axis ij;colormap('gray'),colorbar,

for i=1:nx
   fre1 = 20/(1/dt/2);   %低通19Hz
   [B,A] = butter(8,fre1,'low');
   friraf_bd(1:nt,i) = filtfilt(B,A,riraf_bd(1:nt,i));
end
figure,pcolor(friraf_bd);shading interp; axis ij;colormap('gray'),colorbar,
%%
figure,pcolor(10*shotground+0.001*riraf_bd-0.01*lfils-friraf_bd);shading interp; axis ij;colormap('gray'),colorbar,


