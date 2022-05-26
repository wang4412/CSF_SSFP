% Generate a CSF flow SSFP pattern dictionary
%
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA

%%
clear all;
close all;

acFlow=0:5:30; %mm/s
dcFlow=-5:5; %mm/s

Seq_para.band_dist=22.5; % mm
Seq_para.TR = 6;  % ms
Seq_para.FA = 45; % degree
Seq_para.nvox = 450;  

OPTIONS.tolerance=0.05;
OPTIONS.quiet = 0;
% OPTIONS.nCycle = 20;

Flow_para.T1=4000;  
Flow_para.T2=2000;
Flow_para.vox_length=0.1; % mm

period = 1000; %ms
timeRatio=1; %of out and in flow, >1  

%%
nvox = Seq_para.nvox;
velocity=zeros(length(acFlow),length(timeRatio),length(dcFlow),period);
finalFlow = zeros(length(acFlow),length(timeRatio),length(dcFlow),...
     nvox,floor(period/Seq_para.TR));
finalMag = zeros(3,length(acFlow),length(timeRatio),length(dcFlow),...
     nvox,floor(period/Seq_para.TR));
iters=zeros(length(acFlow),length(timeRatio),length(dcFlow));
rel_diffs=zeros(length(acFlow),length(timeRatio),length(dcFlow));

ElapsedTime = zeros(length(acFlow),length(timeRatio),length(dcFlow));
tic
for i=1:length(acFlow)
    for k=1:length(timeRatio)
        for m=1:length(dcFlow)
            
            fprintf('********************\n');
            fprintf('Processing i,m:%d,%d\n',i,m);

            ParaSet.peakIn=acFlow(i);
            ParaSet.peakOut=acFlow(i);
            ParaSet.timeRatio=timeRatio(k);
            ParaSet.dcFlow=dcFlow(m);
            ParaSet.period=period;

            Flow_para.velocity =FlowVelocity(ParaSet);
            velocity(i,k,m,:)=Flow_para.velocity;

            [flowPattern,rel_diff,magnetizn,~]=...
                SSFP_Flow_Simu(Flow_para,Seq_para,OPTIONS);

            finalFlow(i,k,m,:,:)=flowPattern(:,:,end);% only take the final steady state signal

            iters(i,k,m)=size(flowPattern,3);
            rel_diffs(i,k,m) = rel_diff(end);
            ElapsedTime(i,k,m)=toc;
        end
    end
end

%% DownSampling

tshift_step=8;
timePoint_card = [3 4 5];
yshift=5; % how many pixels to extend
yextract = round((225-13*yshift):13:(225+13*yshift));
lookupFlow = permute(abs(finalFlow),[4 5 1 3 2]);
% imshow3(lookupFlow(125:350,:,:,:),[-0.1 0.6],[length(dcFlow) length(acFlow)]);colorbar;

Dic_down = lookupFlow(yextract,:,:,:);
[dicLength, dicTime,ac_num, dc_num] = size(Dic_down);
tshift_amount = floor(size(Dic_down,2)/tshift_step);
Dic_2down=nan*zeros(dicLength,max(timePoint_card),ac_num, dc_num, tshift_step, length(timePoint_card));
for k=1:length(timePoint_card)
    for i=1:tshift_step
        extract_time = round(linspace(1,size(lookupFlow,2),timePoint_card(k)+1));     
        Dic_2down(:,1:timePoint_card(k),:,:,i,k) = Dic_down(:,extract_time(1:end-1),:,:);
        Dic_down = circshift(Dic_down,-tshift_amount,2);
    end
end
% imshow3(Dic_2down(:,:,:,:,1,end),[-0.1 0.6],[11 7]);colorbar;

%% Vectorize and normalize dictionray entries
dicSize = size(Dic_2down);% dicLength,dicTime,ac_num,dc_num,tshift,tlength
dic_norm = reshape(Dic_2down,prod(size(Dic_2down,[1 2])),[]);
dic_norm = dic_norm-mean(dic_norm,'omitnan');
dic_norm(isnan(dic_norm))=0;
dic_norm = dic_norm./vecnorm(dic_norm);

save dicNorm_dist22p5_y11 Dic_2down dic_norm dicSize tshift_step timePoint_card acFlow dcFlow


