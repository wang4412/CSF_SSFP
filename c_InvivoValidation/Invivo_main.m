% In vivo quantification for CSF dynamics
%
% Created by Yicun Wang (yicun.wang@nih.gov)
% AMRI, LFMI, NINDS, National Institutes of Health, Bethesda, MD, USA
%

clear all;
close all;

load test_data
figure;
imshow(reshape(permute(CSF_patterns,[1 3 2]),[],size(CSF_patterns,2)),[])

%% Check Phyio data
[pks_card,locs_card] = ...
    findpeaks(cardiac,'MinPeakDistance',12);%,'MinPeakHeight',0.8);

figure;
l1=subplot(211);
plot(timeVec,cardiac);
hold on;
plot(timeVec(locs_card),pks_card,'ro');
xlabel 'Time(s)'
l2=subplot(212);plot(timeVec,respiration);
linkaxes([l1,l2],'x');
xlim([0 445])
locs_card_TR = round(timeVec(locs_card)/(TR_ms*1e-3));

diff_locs_card_TR = diff(locs_card_TR);
diff_locs_card_s = diff(timeVec(locs_card));

%% Dictionary matching
load ../a_Dictionary/dicNorm_dist22p5_y11_stepp2.mat

max_tPoints = max(timePoint_card);
[y_extent, t_extent, roi_num] = size(CSF_patterns);
patches = nan*zeros(y_extent,max_tPoints,length(locs_card_TR),roi_num);
tSteps_index = zeros(length(locs_card_TR),1); %Same assignment for all rois

for i=2:length(locs_card_TR)-1    
    if diff_locs_card_s(i)<mean([0.75,1])
        card_range = (locs_card_TR(i)):locs_card_TR(i)+2;
        tSteps_index(i)=3;
        patches(:,1:3,i,:) = CSF_patterns(:,card_range,:);        
    elseif diff_locs_card_s(i)<mean([1,1.25])
        card_range = (locs_card_TR(i)):locs_card_TR(i)+3;
        tSteps_index(i)=4;
        patches(:,1:4,i,:) = CSF_patterns(:,card_range,:);
    else
        card_range = (locs_card_TR(i)):locs_card_TR(i)+4;
        tSteps_index(i)=5;
        patches(:,:,i,:) = CSF_patterns(:,card_range,:);
    end
end

patches_vec = reshape(patches,prod(size(patches,1,2)),prod(size(patches,3,4)));
patches_vec = patches_vec-mean(patches_vec,'omitnan');
patches_vec(isnan(patches_vec))=0;
patches_vec = patches_vec./vecnorm(patches_vec);

corr_Mat = patches_vec'*dic_norm;
corr_Mat=reshape(corr_Mat,[length(locs_card_TR) roi_num dicSize(3:end)]);
[corr_Max idx] = max(corr_Mat,[],3:6,'linear');
[~,~,idxAc,idxDc,idxTshift,idxTsteps] = ind2sub(size(corr_Mat),idx);
ac_results = acFlow(idxAc); 
dc_results = dcFlow(idxDc);

idx_dic = sub2ind(size(Dic_2down,[3,4,5,6]),idxAc,idxDc,idxTshift,idxTsteps);
dic_data = Dic_2down(:,:,idx_dic);
dic_data = reshape(dic_data,[size(dic_data,[1 2]),length(locs_card_TR),roi_num]);

ac_results(corr_Max<0.4)=nan;
dc_results(corr_Max<0.4)=nan;

%% Evaluation

time_range = [0 timeVec(end)];
flow_range = [-3 3];

check_loc = 1;
f1=figure;
c=get(gca,'colororder');

subplot(511);hold on;
l1=plot(timeVec,respiration*8);l1.Color(4)=0.8;
plot(timeVec(locs_card),ac_results(:,check_loc))
xlim(time_range)
ylabel 'AC in mm/s'

subplot(512);hold on;
l2=plot(timeVec,respiration*2);l2.Color(4)=0.8;
plot(timeVec(locs_card),dc_results(:,check_loc)-nanmean(dc_results(:,check_loc)))
xlim(time_range);ylim(flow_range)
plot([time_range(1) time_range(end)],[0 0],'k--');
ylabel 'DC in mm/s'

subplot(513);hold on;
l3=plot(timeVec(locs_card),idxTshift(:,check_loc),'.');
l3.Color=c(2,:);%Orange
xlim(time_range);ylim([1 8])
ylabel 'tShift'

subplot(514);hold on;
l3=plot(timeVec(locs_card),idxTsteps(:,check_loc)+2,'o');
l3.Color=c(2,:);%Orange
xlim(time_range);
ylabel 'card t steps'; 

subplot(515);hold on;
l3=plot(timeVec(locs_card),corr_Max(:,check_loc));
l3.Color=c(2,:);%Orange
xlim(time_range)
ylim([0.4 1.2])
ylabel 'corrCoef'; 
xlabel Time(s)


