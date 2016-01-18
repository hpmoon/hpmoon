function [Feat Class]=gentofen_class_jac_g(seg,lev,det_an,Feat);
global class_essex;
global data_essex;
%How to load the data. It has to be executed in the same directory as the
%datafiles 
% seg (segment) goes from 1 to 20
% lev (level) goes from 1 to 7
% electrodes are 15
electrodes=size(data_essex,3);
% det_an=0 means Detail and det_an=1 means Analysis

XX(1:size(data_essex,5),1:electrodes)=0;

for kk = 1 : size(data_essex,5)
    for jj = 1 : electrodes
       XX(kk,jj)=data_essex(seg,lev,jj,det_an+1,kk);
    end   
end

Feat = [Feat XX];
Class = class_essex;

end