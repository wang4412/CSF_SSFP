function [flowPattern rel_diff magnetizn_pattern zposition] ...
    = SSFP_Flow_Simu(Flow_para,Seq_para,OPTIONS)
% Simulate SSFP CSF flow moving in z linear gradient (1ms steps)
%   Output time resolution is TR
%   |j=1                    |j=2                    |j=3
%              TR1                      TR2                    TR3 
%               RF                      RF            
%   | 1 | 2 | 3 | 4 | 5 | 6 | 1 | 2 | 3 | 4 | 5 | 6 | 1 | 2 | 3 |...
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
%               * velocity in mm/s. 1000x1 vector for 1 second with 1ms resolution
%               * chemicalShiftHz
%
% * Seq_para:   * TR in ms. Default 6 ms
%               * FA in degrees. Default 45
%               * band_dist in mm. Default 15
%               * nvox. Total voxel number along the aqueduct. Default 3*band_dist/vox_length
%
% * OPTIONS:    * nCycle. maximum 1s cycles simulated. Default 10
%               * tolerance. RMS change from last cycle to stop. Default
%               0.05
%               * quiet. suppress output TR No.
%   
% ---------------------------
% OUTPUTS:
% ---------------------------
% * flowPattern: zposition corrected Mxy flow pattern. (length,cycle number,TR in each cycle)
% * magnetizn:   zposition corrected simulation results.
%                 ([Mx,My,Mz],length,cycle number,TR in each cycle)
% * rel_dif:     improvement vector
% * zposition:   (length,cycle number,TR in each cycle)
% 
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%

%% Parse input and initialization 
% flow
T1 = 3000; T2 = 1500; vox_length = 0.1; 
velocity = Flow_para.velocity;
chemicalShiftHz = 0;

if isfield(Flow_para,'T1')
    T1 = Flow_para.T1;
end
if isfield(Flow_para,'T2')
    T2 = Flow_para.T2;
end  
if isfield(Flow_para,'vox_length')
    vox_length = Flow_para.vox_length;
end   
if isfield(Flow_para,'chemicalShiftHz')
    chemicalShiftHz = Flow_para.chemicalShiftHz;
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
nTR_inCycle = floor(length(velocity)/TR);

% OPTIONS
nCycle = 10;
if isfield(OPTIONS,'nCycle') % max simulated cycles
    nCycle = OPTIONS.nCycle;
end
if isfield(OPTIONS,'tolerance')
    tolerance = OPTIONS.tolerance;
end
if isfield(OPTIONS,'quiet')
    quiet = OPTIONS.quiet;
end

% Initializtion Calcs
gradz = 1e3/TR/42.6/band_dist; %mT/m

zposition = zeros(nvox,nCycle,nTR_inCycle);% (length,cycle number,TR in each cycle)
zposition_temp = (-(1:nvox)*vox_length)';

magnetizn = zeros(3,nvox,nCycle,nTR_inCycle);% ([Mx,My,Mz],length,cycle number,TR in each cycle)
M_temp = zeros(3,nvox);
M_temp(3,:)=1;
magnetizn_pattern = zeros(3,nvox,nCycle,nTR_inCycle);% ([Mx,My,Mz],length,cycle number,TR counter in each cycle)

flowPattern = zeros(nvox,nCycle,nTR_inCycle);% (length,cycle number,TR in each cycle)
rel_diff = zeros(nCycle,1);

excitationRotMat =  Bloch_FA_RotMat(FA,0,90);% B1 along x axis
RelaxMat1 = diag([exp(-1/T2) exp(-1/T2) exp(-1/T1)]); % Mx My Mz for each 1 ms!
RelaxMat2 = [0 0 1-exp(-1/T1)]';

%% Simulation

for i=1:nCycle   % over cardiac cycles
    if ~quiet
        fprintf('Processing Cycle %d\n',i)
        fprintf('TR');
    end
    for j=1:nTR_inCycle % indexing 6ms cycles within 1s
        if ~quiet
            fprintf('%d..',j);
        end
        
        for k = 1:TR % 1 ms time steps in TR 
            wot = squeeze(zposition_temp*gradz*42.6*2*pi*1e-3...
                +chemicalShiftHz*2*pi*1e-3); %rad  
            parfor m = 1:nvox % voxels along the length
                % Phase change due to flow               
                M_temp(:,m) = zRotMat(wot(m))*M_temp(:,m);
                % T1 T2 Relaxation
                M_temp(:,m) = RelaxMat1*M_temp(:,m)+RelaxMat2;
                % RF
                if k==round(TR/2)
                    M_temp(:,m) = excitationRotMat*M_temp(:,m);
                end
            end 
            velocity_temp = velocity((j-1)*TR+k);           
            zposition_temp = zposition_temp+velocity_temp*1e-3;  % mm               
        end
            
        magnetizn(:,:,i,j) = M_temp;
        zposition(:,i,j) = zposition_temp;  
        
        Shift_pix = round((zposition_temp(1)-zposition(1,1,1))./vox_length);
        magnetizn_pattern(:,:,i,j) = circshift(sz(magnetizn(:,:,i,j)),-Shift_pix,2); 
    end    
    
    Mxy = sz(magnetizn_pattern(1,:,i,:)+1i*magnetizn_pattern(2,:,i,:));
    flowPattern(:,i,:) = abs(Mxy);
    
    if i==1
        comp = zeros(size(flowPattern(:,i,:)));
    else
        comp = flowPattern(:,i-1,:);
    end
    rel_diff(i) = norm(reshape(flowPattern(:,i,:)-comp,1,[]))...
        /norm(reshape(flowPattern(:,i,:),1,[])); 
    
    if ~quiet
        fprintf('\n');
        fprintf('RelativeDiff %2.2f%%\n',rel_diff(i)*100);
        fprintf('\n');
    end
       
    if exist('tolerance','var')
        if rel_diff(i)<tolerance % Reaching tolerence then create iters
            iters = i;
            break;
        end
    end
end

flowPattern = permute(flowPattern,[1 3 2]);
magnetizn_pattern = permute(magnetizn_pattern,[1 2 4 3]);

if exist('iters','var') % if exitting early, remove redundant space
    flowPattern(:,:,iters+1:end)=[];
    rel_diff(iters+1:end)=[];
    magnetizn_pattern(:,:,:,iters+1:end)=[];
end


