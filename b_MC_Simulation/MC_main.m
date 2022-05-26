% Monte-Carlo Simulation
% Generate 3 freq superimposed flow velocity with random amplitude
% And Decode flows using correlation
%
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%

clear all;
close all;

trial_no = 5;
%%
cardPeriod_num = 300;
respPeriod_num = round(cardPeriod_num/6);
autoPeriod_num = round(cardPeriod_num/20);

Periods = cell(trial_no,3);
Amps = cell(trial_no,3);
Velocities = cell(trial_no,4);

for j=1:trial_no
    cardPeriod =  1000+round(200*randn(cardPeriod_num,1));%ms
    respPeriod =  6000+round(1000*randn(respPeriod_num,1));
    autoPeriod = 20000+round(5000*randn(autoPeriod_num,1));
    %                               std  mean
    cardAmp = randn(cardPeriod_num,1)*3  +6; cardAmp(cardAmp<0)=0;
    respAmp = randn(respPeriod_num,1)*1  +1;  respAmp(respAmp<0)=0;
    autoAmp = randn(autoPeriod_num,1)*1  +3;

    v1=[];v2=[];v3=[];
    for i=1:cardPeriod_num
        v1=[v1;cardAmp(i)*sin((1:cardPeriod(i))'/cardPeriod(i)*2*pi)];
    end
    for i=1:respPeriod_num
        v2=[v2;respAmp(i)*sin((1:respPeriod(i))'/respPeriod(i)*2*pi)];
    end
    for i=1:autoPeriod_num
        v3=[v3;autoAmp(i)*sin((1:autoPeriod(i))'/autoPeriod(i)*2*pi)];
    end

    minLength = min([length(v1),length(v2),length(v3)]);
    v1 = v1(1:minLength);v2 = v2(1:minLength);v3 = v3(1:minLength);
    v = v1+v2+v3;
    
    Periods(j,:) = [{cardPeriod} {respPeriod} {autoPeriod}];
    Amps(j,:) =    [{cardAmp} {respAmp} {autoAmp}];
    Velocities(j,:) = [{v1} {v2} {v3} {v}];
    
    j
end

Seq_para.band_dist=22.5;
Seq_para.TR = 6;
Seq_para.FA = 45;
Seq_para.nvox = 450;
OPTIONS.quiet = 0;

Flow_para.T1=4000;  
Flow_para.T2=2000;
Flow_para.vox_length=0.1;

yshift=5;
y_reso = 1.33/Flow_para.vox_length;
yextract = round((225-y_reso*yshift):y_reso:(225+y_reso*yshift));

FlowPattern = cell(trial_no,1);
tic
for j=1:trial_no
    Flow_para.velocity = Velocities{j,4};
    temp = SSFP_Flow_Simu_Aperiodic(Flow_para,Seq_para,OPTIONS);
    FlowPattern{j} = temp(yextract,123:246:end);
    j
end
T_elapsed = toc;

i=1;
temp = FlowPattern{i};
figure;
plot(Velocities{i,2}+Velocities{i,3})
figure;
imshow(temp,[])

save('velocity_3freq_bulk_500Runs','Periods','Amps','Velocities',...
    'FlowPattern','T_elapsed','yshift','Seq_para','Flow_para','trial_no');

%%

load ../a_Dictionary/dicNorm_dist22p5_y11

yshift=5;
y_extent = 2*yshift+1;
TR_ms = 246;

dcFlow_idx = 1:length(dcFlow);
acFlow_idx = 1:length(acFlow);
dic_norm = reshape(dic_norm,dicSize);
dic_norm = dic_norm(:,:,acFlow_idx,dcFlow_idx,:,:);
dic_norm = reshape(dic_norm,prod(dicSize(1:2)),[]);
dicSize(3) = length(acFlow_idx);
dicSize(4) = length(dcFlow_idx);

Flow_results = cell(trial_no,6);
patches_all = cell(trial_no);
tSteps_index_all = cell(trial_no);

for j = 1:trial_no
    cardPeriod = Periods{j,1};
    v = Velocities{j,4};
    flowPattern = FlowPattern{j};
    
    locs_card = round(cumsum(cardPeriod)-3/4*cardPeriod);
    locs_card(locs_card>length(v))=[];    
    diff_locs_card_s = cardPeriod/1000;
    locs_card_TR = round(locs_card/TR_ms);
    csf_roi = flowPattern;
    
    max_tPoints = max(timePoint_card);
    patches = nan*zeros(y_extent,max_tPoints,length(locs_card_TR));
    tSteps_index = zeros(length(locs_card_TR),1); %Same assignment for all rois

    for i=2:length(locs_card_TR)-1    
        if diff_locs_card_s(i)<mean([0.75,1])
            range = (locs_card_TR(i)):locs_card_TR(i)+2;
            tSteps_index(i)=3;
            patches(:,1:3,i) = csf_roi(:,range,:,:);        
        elseif diff_locs_card_s(i)<mean([1,1.25])
            range = (locs_card_TR(i)):locs_card_TR(i)+3;
            tSteps_index(i)=4;
            patches(:,1:4,i) = csf_roi(:,range,:,:);
        else
            range = (locs_card_TR(i)):locs_card_TR(i)+4;
            tSteps_index(i)=5;
            patches(:,:,i) = csf_roi(:,range,:,:);
        end
    end

    patches_vec = combine_dim(patches,[1 2]);
    patches_vec = patches_vec-mean(patches_vec,'omitnan');
    patches_vec(isnan(patches_vec))=0;
    patches_vec = patches_vec./vecnorm(patches_vec);

    corr_Mat = patches_vec'*dic_norm;
    corr_Mat=reshape(corr_Mat,[length(locs_card_TR) dicSize(3:end)]);

    [corr_Max idx] = max(corr_Mat,[],2:5,'linear');
    [~,idxAc,idxDc,idxTshift,idxTsteps] = ind2sub(size(corr_Mat),idx);
    
    ac_results = acFlow(idxAc); dc_results = dcFlow(idxDc);  
    Flow_results(j,:) = [{ac_results} {dc_results} {idxTshift} {idxTsteps}...
        {corr_Max},{locs_card}];
    patches_all{j} = patches;
    tSteps_index_all{j} = tSteps_index;
end

save('velocity_DetResults',...
    'Flow_results','patches_all','tSteps_index_all');%'SNR','noise_std',

%%
j=3;

ac_results= Flow_results{j,1};
dc_results= Flow_results{j,2};
idxTshift= Flow_results{j,3};
idxTsteps = Flow_results{j,4};
corr_Max = Flow_results{j,5};
locs_card = Flow_results{j,6};

flowPattern = FlowPattern{j};
cardAmp = Amps{j,1};
v1 = Velocities{j,1};
v2 = Velocities{j,2};
v3 = Velocities{j,3};
v = Velocities{j,4};
tSteps_index = tSteps_index_all{j};
idx_err=find(abs(dc_results'-(v2(locs_card)+v3(locs_card)))>5);

figure('position',[0    0.2633    0.4359    0.5342]);
subplot(411);hold off
p1=plot(locs_card/1000,cardAmp(1:length(locs_card)));hold on
p2=plot(locs_card/1000,ac_results);
p1.Color(4)=0.8;p2.Color(4)=0.8;
xlim([0 300]);box off;
ylim([0 20]);
yticks([0 20])
ax1=get(gca);
ax1.XAxis.Visible = 'off';

subplot(412);hold off
p1=plot(locs_card/1000,v2(locs_card)+v3(locs_card));hold on;
p2=plot(locs_card/1000,dc_results);
ylim([-10 10]);
xlim([0 300]);box off;
p1.Color(4)=0.8;p2.Color(4)=0.8;
ax1=get(gca);
ax1.XAxis.Visible = 'off';

subplot(413);
p1=plot(locs_card/1000,tSteps_index,'o');hold on
p2=plot(locs_card/1000,idxTsteps+2,'.','markersize',10);
p1.Color(4)=0.8;p2.Color(4)=0.8;
ylim([2.5 5.5])
xlim([0 300]);box off;
ax1=get(gca);
ax1.XAxis.Visible = 'off';

subplot(414);
plot(locs_card/1000,corr_Max,'k');
ylim([0.6 1.1]);
xlim([0 300]);box off;
xticks(0:50:300)

