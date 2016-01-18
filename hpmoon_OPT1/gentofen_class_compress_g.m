function [Feat Class]=gentofen_class_compress_g(seg,lev,elect,det_an,Feat);
global class_essex;
global data_essex;
%It has to be executed in the same directory as the datafiles 
% seg (segment) goes from 1 to 20
% lev (level) goes from 1 to 7
% elect (electrode) goes from 1 to 15
% det_an=0 means Detail and det_an=1 means Analysis

for kk = 1 : size(data_essex,5)
       XX(kk,1)=data_essex(seg,lev,elect,det_an+1,kk);
end


Feat = [Feat XX];
Class = class_essex;

end