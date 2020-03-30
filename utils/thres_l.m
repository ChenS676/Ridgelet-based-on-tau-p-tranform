%% 
function [bd_reccoef,r,c] = thres_l(s,mode,scales)
[r,c] = size(s);    %在radon域滤波选择这个

% part_l(:,i) = reshape(s,nt*45,1); 
switch mode
    case 'gauss'
        gr_reccoef = rifilter2(s,r,c);

    case 'linear'
        fre = 30;       
        scale = 1/fre;
        bd_reccoef = (scales < scale).*s;

    case 'parabolic'      %实际数据gmdground
    filts = mwthighp(part_l,r,c,nx)
end 