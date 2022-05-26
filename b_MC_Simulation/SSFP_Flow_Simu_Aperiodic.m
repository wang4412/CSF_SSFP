function [flowPattern magnetizn_pattern] = ...
    SSFP_Flow_Simu_Aperiodic(Flow_para,Seq_para,OPTIONS)
% Simulate aperiodic CSF flow 
% ---------------------------
% TEST:
% ---------------------------
% flow_para.velocity=20*ones(1000,1);
% flowPattern=SSFP_Flow_Simu(flow_para,[],[]);
% figure;imshow(reshape(flowPattern,450,[]),[])
%
% ---------------------------
% INPUTS:
% ---------------------------
% * Flow_para:  * T1 in ms. Default 3000
%               * T2 in ms. Default 1500
%               * vox_length in mm. Single vox length. Default 0.1
%               * velocity in mm/s. vector.
%               
% * Seq_para:   * TR in ms. Default 6 ms
%               * FA in degrees. Default 45
%               * band_dist in mm. Default 15
%               * nvox. Total voxel number along the aqueduct. Default 3*band_dist/vox_length
%
% * OPTIONS:    * quiet. suppress output TR No.
%   
% ---------------------------
% OUTPUTS:
% ---------------------------
% * flowPattern: zposition corrected Mxy flow pattern. (length,cycle number,TR in each cycle)
% 
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA

%% Parse input and initialization 
lineLength=0;

% flow
T1 = 3000; T2 = 1500; vox_length = 0.1; 
velocity = Flow_para.velocity;

if isfield(Flow_para,'T1')
    T1 = Flow_para.T1;
end
if isfield(Flow_para,'T2')
    T2 = Flow_para.T2;
end  
if isfield(Flow_para,'vox_length')
    vox_length = Flow_para.vox_length;
end         
    
% Sequence
TR = 6; FA = 45; band_dist = 15;
if isfield(Seq_para,'TR')
    TR = Seq_para.TR;
end
if isfield(Seq_para,'FA')
    FA = Seq_para.FA;
end
if isfield(Seq_para,'band_dist')
    band_dist = Seq_para.band_dist;
end
if isfield(Seq_para,'nvox')
    nvox = Seq_para.nvox;
else
    nvox = round(3*band_dist/vox_length);
end 
nMS = length(velocity);

% Initializtion Calcs
gradz = 1e3/TR/42.6/band_dist; %mT/m

zposition = zeros(nvox,nMS);% (length,TR in each cycle)
zposition(:,1) = -(1:nvox)*vox_length;

magnetizn = zeros(3,nvox,nMS);% ([Mx,My,Mz],length,TR in each cycle)
magnetizn(3,:,1) = 1;
M_temp = zeros(3,nvox);
M_temp(3,:)=1;

magnetizn_pattern = zeros(3,nvox,nMS);% ([Mx,My,Mz],length,TR in each cycle)
magnetizn_pattern(3,:,1) = 1;

excitationRotMat =  Bloch_FA_RotMat(FA,0,90);% B1 along x axis
RelaxMat1 = diag([exp(-1/T2) exp(-1/T2) exp(-1/T1)]); % Mx My Mz
RelaxMat2 = [0 0 1-exp(-1/T1)]';

%% Simulation

for i=2:nMS      
    zposition(:,i)= zposition(:,i-1)+velocity(i)*1e-3;  % mm
    wot = zposition(:,i)*gradz*42.6*2*pi*1e-3; % rad
    
    parfor m = 1:nvox % voxels along the length
        % Phase change due to flow               
        M_temp(:,m) = zRotMat(wot(m))*M_temp(:,m);
        % T1 T2 Relaxation
        M_temp(:,m) = RelaxMat1*M_temp(:,m)+RelaxMat2;
        % RF
        if mod(i,TR)==round(TR/2)
            M_temp(:,m) = excitationRotMat*M_temp(:,m);
        end
    end 
    
    magnetizn(:,:,i) = M_temp;
    
    if OPTIONS.quiet~=1 && mod(i,500)==0
        fprintf(repmat('\b', 1, lineLength));
        lineLength=fprintf('%d%%\n',round(100*i/nMS));
    end
    
end   

Shift_pix = squeeze(round((zposition(1,:)-repmat(zposition(1,1),...
    [1 nMS]))./vox_length));
for j = 1:nMS
    magnetizn_pattern(:,:,j) = circshift(sz(magnetizn(:,:,j)),-round(Shift_pix(j)),2);
end   

flowPattern = sz(abs(magnetizn_pattern(1,:,:)+1i*magnetizn_pattern(2,:,:)));


