%% 
function [bd_reccoef,r,c] = thres_l(s,mode,scales)
[r,c] = size(s);    %��radon���˲�ѡ�����

% part_l(:,i) = reshape(s,nt*45,1); 
switch mode
    case 'gauss'
        gr_reccoef = rifilter2(s,r,c);

    case 'linear'
        fre = 30;       
        scale = 1/fre;
        bd_reccoef = (scales < scale).*s;

    case 'parabolic'      %ʵ������gmdground
    filts = mwthighp(part_l,r,c,nx)
end 