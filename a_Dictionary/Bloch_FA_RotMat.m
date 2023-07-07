function mat = Bloch_FA_RotMat(alpha,theta, ksi)
%Genertaing Rotational matrix for Bloch simulation
%   Angles in degree
% 
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%
%%

    alpha = alpha/180*pi;
    theta = theta/180*pi;
    ksi   = ksi/180*pi;
    
    Rz_theta = [cos(theta)      sin(theta)      0;...
                -sin(theta)     cos(theta)      0;...
                0               0               1];
       
    Ry_ksi   = [cos(ksi)        0               -sin(ksi);...
                0               1               0       ;...
                sin(ksi)        0               cos(ksi)];    
            
    Rz_alpha = [cos(alpha)      sin(alpha)      0;...
                -sin(alpha)     cos(alpha)      0;...
                0               0               1];            
    
    mat = inv(Rz_theta)*inv(Ry_ksi)*Rz_alpha*Ry_ksi*Rz_theta;

end

