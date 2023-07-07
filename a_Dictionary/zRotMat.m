function mat = zRotMat(wot)
% Generate Rotational Matrix 3x3xlength(wot)
%
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA

%%
%     mat = [cos(wot) sin(wot)  0;
%           -sin(wot) cos(wot)  0;
%                  0         0  1];

    mat = zeros(3,3,length(wot));
    mat(1,1,:) = cos(wot);
    mat(2,1,:) = -sin(wot);
    mat(1,2,:) = sin(wot);
    mat(2,2,:) = cos(wot);
    mat(3,3,:) = 1;
    
end