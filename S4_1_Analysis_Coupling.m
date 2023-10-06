clc;clear;
couplingpath='D:/Work/hcp_development/results/coupling-multimodal/glm/210/';
savepath='D:\Work\hcp_development\results\coupling-multimodal\FigFiles\';
bna2yeo=load('bna2yeo.mat').bna2yeo;
path='F:/fmri2win';
unrelated=1;
atlas=210;
nums=500;
euc=0;
if euc==0
    dim=4;
else
    dim=5;
end
%% remove fmri outlier subjects
fmripath='D:/Work/hcp_development/results/network/fmri-fisher/246';
subs=dir(fmripath);
subs(1:2)=[];
num=1;
for i=1:length(subs)
    id=split(subs(i).name,'.');
    temp1=load([path,'/',id{1},'_V1_MR/censored_mask.mat']).censored_mask;
    temp2=load([path,'/',id{1},'_V1_MR/Movement_RelativeRMS.mat']).Movement_RelativeRMS;
    temp3=load([path,'/',id{1},'_V1_MR/Movement_RelativeRMS_mean.mat']).Movement_RelativeRMS_mean;
    fram_per=ceil(size(temp1,2)/size(temp3,2));
%     for j=1:size(temp3,2)
%         if temp1(1+fram_per*(j-1):min(fram_per*j))
%     end
    if sum(temp1) > 4*60/0.8
        fmri_id{num,1}=id{1};
        Movement_RelativeRMS(num,1)=sum(temp1);
        Movement_RelativeRMS(num,2)=sum(temp1)/size(temp1,2);
        Movement_RelativeRMS(num,3)=mean(temp2(temp1==1));
        Movement_RelativeRMS(num,4)=std(temp2(temp1==1));
        num=num+1;
    end
end
subs=dir(couplingpath);subs(1:2)=[];
num=1;
for i=1:length(subs)
    for j=1:size(fmri_id,1)
        if strcmp(subs(i).name,[fmri_id{j},'.mat'])
            id_corrected{num,1}=subs(i).name;
            num=num+1;
            break
        end
    end
end
if unrelated==1
    [data,ids,~]=xlsread('D:\Work\hcp_development\adni_babri.xlsx','HCP_D');
    ids(1,:)=[];
    j=1;
    k=1;
    for s=1:length(id_corrected)
        temp=split(id_corrected{s},'.mat');
        for i=j:size(ids,1)
            if strcmp(ids{i,1},temp{1}) &&  data(i,4)==0
                outindex(k)=s;
                k=k+1;
                j=i;
                break;
            end
        end
    end
    id_corrected(outindex,:)=[];
end
save('D:\Work\hcp_development\results\coupling-multimodal\id_corrected.mat','id_corrected');

R=zeros(atlas,dim*2+6,length(id_corrected));
R_null=zeros(atlas,dim*2+6,length(id_corrected),nums);

for s=1:length(id_corrected)
    % w; b; r2; F value; p value; wuchafangcha; adjustedr2; PFM
    temp=load([couplingpath,id_corrected{s}]).resutls;
    R(:,:,s)=temp(:,:,1);
    R_null(:,:,s,:)=temp(:,:,2:nums+1);
end
%% global adjustedr2 cross 246 regions (8 yeo network)
adjusted_r2=squeeze(R(:,dim+6,:));
save([savepath,'adjusted_r2.mat'],'adjusted_r2');
adjusted_r2_null=squeeze(R_null(:,dim+6,:,1:nums));
mask = ~isnan(adjusted_r2) & ~isinf(adjusted_r2);
adjusted_r2(~mask)=nan;
adjusted_r2_mean=nanmean(adjusted_r2,2);
adjusted_r2_null_mean=squeeze(nanmean(adjusted_r2_null,2));
adjusted_r2_mean_p=sum(adjusted_r2_mean<adjusted_r2_null_mean,2)/nums;
adjusted_r2_mean_p_fdr=mafdr(adjusted_r2_mean_p,'BHFDR', true);
adjusted_r2_mean_p_bonf=adjusted_r2_mean_p*atlas;
adjusted_r2_global = [adjusted_r2_mean adjusted_r2_mean_p adjusted_r2_mean_p_fdr adjusted_r2_mean_p_bonf]; % mean adjusted_r2; p
save([savepath,'adjusted_r2_global.mat'],'adjusted_r2_global');
adjusted_r2_global_yeo=zeros(max(bna2yeo(:,1)),2);
temp=0;
adjusted_r2_global_yeo_anova=[0 0];
for i=1:max(bna2yeo(:,1))
    adjusted_r2_global_yeo(i,1)=mean(adjusted_r2_mean(bna2yeo(:,1)==i));
    adjusted_r2_global_yeo(i,2)=sum(mean(adjusted_r2_mean(bna2yeo(:,1)==i))<mean(adjusted_r2_null_mean(bna2yeo(:,1)==i,:)))/nums;
    adjusted_r2_global_yeo_draw{i}=adjusted_r2_mean(bna2yeo(:,1)==i);
    temp=temp+sum(adjusted_r2_mean(bna2yeo(:,1)==i));
    adjusted_r2_global_yeo_anova=[adjusted_r2_global_yeo_anova;[adjusted_r2_global_yeo_draw{i} ones(size(adjusted_r2_global_yeo_draw{i}))*i]];
end
adjusted_r2_global_yeo_anova(1,:)=[];
[~,tb1,stats]=anova1(adjusted_r2_global_yeo_anova(:,1),adjusted_r2_global_yeo_anova(:,2),'off');
[post_hoc,~,~,~]=multcompare(stats);
temp/atlas;

[~,tb1,stats1]=kruskalwallis(adjusted_r2_global_yeo_anova(:,1),adjusted_r2_global_yeo_anova(:,2),'off');
[post_hoc1,~,~,~]=multcompare(stats1);
% for i=1:size(post_hoc,1)
%     post_hoc_p(post_hoc(i,1),post_hoc(i,2))=post_hoc(i,6);
%     post_hoc_p(post_hoc(i,2),post_hoc(i,1))=post_hoc(i,6);
% end
adjusted_r2_global_yeo(:,3)=mafdr(adjusted_r2_global_yeo(:,2),'BHFDR', true);
adjusted_r2_global_yeo(:,4)=adjusted_r2_global_yeo(:,2)*atlas;
save([savepath,'adjusted_r2_global_yeo.mat'],'adjusted_r2_global_yeo'); % mean adjusted_r2; p

%% significant number fraction across 246 regions (8 yeo network)
% p=squeeze(R(:,dim+4,:));
p=sum(adjusted_r2<adjusted_r2_null,3)/nums;% permutation using null models
sigfrac(:,1)=sum(p<0.05,2)./size(p,2);
p_fdr=ones(size(p))*999;
temp=mafdr(p(mask),'BHFDR', true);
p_fdr(mask)=temp;
sigfrac(:,2)=sum(p_fdr<0.05,2)./size(p,2);
p_bonf=ones(size(p))*999;
temp=p(mask)*size(p(mask),1);
p_bonf(mask)=temp;
sigfrac(:,3)=sum(p_bonf<0.05,2)./size(p,2);
save([savepath,'sigfrac.mat'],'sigfrac');

adjusted_r2_yeo=zeros(max(bna2yeo(:,1)),size(adjusted_r2,2));
adjusted_r2_null_yeo=zeros(max(bna2yeo(:,1)),size(adjusted_r2,2),nums);
for i=1:max(bna2yeo(:,1))
    adjusted_r2_yeo(i,:)=mean(adjusted_r2(bna2yeo(:,1)==i,:));
    adjusted_r2_null_yeo(i,:,:)=mean(adjusted_r2_null(bna2yeo(:,1)==i,:,:));
end
p_yeo=sum(adjusted_r2_yeo<adjusted_r2_null_yeo,3)/nums;% permutation using null models
sigfrac_yeo(:,1)=sum(p_yeo<0.05,2)./sum(p_yeo,2);
sigfrac(:,1)=sum(p<0.05,2)./size(p,2);
p_yeo_fdr=ones(size(p_yeo))*999;
mask_yeo = ~isnan(adjusted_r2_yeo) & ~isinf(adjusted_r2_yeo);
temp=mafdr(p_yeo(mask_yeo),'BHFDR', true);
p_yeo_fdr(mask_yeo)=temp;
sigfrac_yeo(:,2)=sum(p_yeo_fdr<0.05 & p_yeo_fdr~=999,2)./sum(p_yeo_fdr~=999,2);

p_yeo_bonf=ones(size(p_yeo))*999;
temp=p_yeo(mask_yeo)*size(p_yeo(mask_yeo),1);
p_yeo_bonf(mask_yeo)=temp;
sigfrac_yeo(:,3)=sum(p_yeo_bonf<0.05 & p_yeo_bonf~=999,2)./sum(p_yeo_bonf~=999,2);
save([savepath,'sigfrac_yeo.mat'],'sigfrac_yeo');
%% PFM weight of property across 246 regions (8 yeo network)
w_pfm=R(:,end-dim+1:end,:);
w_pfm_null=R_null(:,end-dim+1:end,:,:);
save([savepath,'w_pfm.mat'],'w_pfm');
w_pfm_mean=nanmean(w_pfm,3);
save([savepath,'w_pfm_mean.mat'],'w_pfm_mean');
w_pfm_null_mean=squeeze(nanmean(w_pfm_null,3));
w_pfm_mean_p=sum(abs(w_pfm_mean)<abs(w_pfm_null_mean),3)/nums;
w_pfm_mean_corr=corr(w_pfm_mean,'type','pearson');
save([savepath,'w_pfm_mean_p.mat'],'w_pfm_mean_p');
w_pfm_mean_p_fdr=mafdr(w_pfm_mean_p(:),'BHFDR', true);
w_pfm_mean_p_fdr=reshape(w_pfm_mean_p_fdr,size(w_pfm_mean_p));
save([savepath,'w_pfm_mean_p_fdr.mat'],'w_pfm_mean_p_fdr');
w_pfm_mean_p_bonf=w_pfm_mean_p*size(w_pfm_mean_p,1)*size(w_pfm_mean_p,2);
w_pfm_mean_p_bonf=reshape(w_pfm_mean_p_bonf,size(w_pfm_mean_p));
save([savepath,'w_pfm_mean_p_bonf.mat'],'w_pfm_mean_p_bonf');

% feature prefer to predict different regions
w_pfm_prefer_mask=zeros(size(w_pfm_mean));
temp2=zeros(size(w_pfm,2)*size(w_pfm,3),2);
for i=1:atlas
    temp1=squeeze(w_pfm(i,:,:));
    temp1(end-2,:)=-(temp1(end-2,:));
%     for j=1:size(temp1,1)
%         signs(i,j)=ttest2(temp1(1,:),abs(temp1(1,:)));
%     end
%     temp1(:,signs(i,:)>0.05)=abs(temp1(:,signs(i,:)>0.05));
    temp2(:,1)=[temp1(1,:)';temp1(2,:)';temp1(3,:)';temp1(4,:)'];
    temp2(:,2)=[ones(size(temp1,2),1);ones(size(temp1,2),1)*2;ones(size(temp1,2),1)*3;ones(size(temp1,2),1)*4];
    temp2=abs(temp2);
    [~,w_pfm_prefer_tb1{i,1},w_pfm_prefer_stats{i,1}]=anova1(temp2(:,1),temp2(:,2),'off');
    [w_pfm_prefer_post_hoc{i,1},~,~,~]=multcompare(w_pfm_prefer_stats{i,1});
    w_pfm_prefer_p(i,1)=w_pfm_prefer_tb1{i, 1}{2, 6};
    [~,index]=sort(w_pfm_prefer_stats{i, 1}.means,'descend');
    flag=find(w_pfm_prefer_post_hoc{i, 1}(:,1)==min(index(1),index(2)) & w_pfm_prefer_post_hoc{i, 1}(:,2)==max(index(1),index(2)));
    if w_pfm_prefer_post_hoc{i, 1}(flag,6)<=0.05
        w_pfm_prefer_mask(i,index(1))=1;
    else
        flag1=find(w_pfm_prefer_post_hoc{i, 1}(:,1)==min(index(2),index(3)) & w_pfm_prefer_post_hoc{i, 1}(:,2)==max(index(2),index(3)));
        flag2=find(w_pfm_prefer_post_hoc{i, 1}(:,1)==min(index(1),index(3)) & w_pfm_prefer_post_hoc{i, 1}(:,2)==max(index(1),index(3)));
        if w_pfm_prefer_post_hoc{i, 1}(flag1,6)<=0.05
            w_pfm_prefer_mask(i,index(1))=1;
            w_pfm_prefer_mask(i,index(2))=1;
        else
            if w_pfm_prefer_post_hoc{i, 1}(flag2,6)<=0.05
                w_pfm_prefer_mask(i,index(1))=1;
                w_pfm_prefer_mask(i,index(2))=1;
            else
                w_pfm_prefer_mask(i,index(1))=1;
                w_pfm_prefer_mask(i,index(2))=1;
                w_pfm_prefer_mask(i,index(3))=1;
            end
        end
    end
end
w_pfm_prefer_p=mafdr(w_pfm_prefer_p,'BHFDR', true);
w_pfm_prefer_mask(w_pfm_prefer_p>0.05,:)=0;
w_pfm_prefer_mask(5,4)=0;w_pfm_prefer_mask(6,4)=0;w_pfm_prefer_mask(147,:)=[1,1,1,0];%¾À´í
save([savepath,'w_pfm_prefer_mask.mat'],'w_pfm_prefer_mask');
w_pfm_prefer=w_pfm_mean.*w_pfm_prefer_mask;
w_pfm_prefer(w_pfm_prefer_p>0.05,:)=0;
save([savepath,'w_pfm_prefer.mat'],'w_pfm_prefer');

for i=1:max(bna2yeo(:,1))
    w_pfm_prefer_yeo(i,:)=sum(w_pfm_prefer_mask(bna2yeo(:,1)==i,:))/sum(bna2yeo(:,1)>0);
end

[~,index]= max(abs(w_pfm_mean),[],2);
[~,index_null]= max(abs(w_pfm_null_mean),[],2);
index_null=squeeze(index_null);
if euc==0
    index=index+1;
    index_null=index_null+1;
end
ind_stat=tabulate(index);
ind_null_stat=zeros(size(ind_stat,1),size(ind_stat,2),nums);
for i=1:nums
    ind_null_stat(:,:,i)=tabulate(index_null(:,i));
end
ind_stat(:,4)=sum(ind_stat(:,2)<squeeze(ind_null_stat(:,2,:)),2)/nums;
save([savepath,'ind_stat.mat'],'ind_stat');
%% coupling alteration with age
[data,ids,~]=xlsread('D:\Work\hcp_development\adni_babri.xlsx','HCP_D');
% [data,ids,~]=xlsread('D:\Work\hcp_development\hcd_cognitivescore_results.xlsx','Results');
ids(1,:)=[];
demography=zeros(length(id_corrected),size(data,2));
j=1;
for s=1:length(id_corrected)
    temp=split(id_corrected{s},'.mat');
    for i=j:size(ids,1)
        if strcmp(ids{i,1},temp{1})
            demography(s,:)=data(i,:);
            j=i;
            tbv(s,:)=load([path,'/',ids{i,1},'_V1_MR/tbv.mat']).tbv;
%             gs_surf(s,:)=load([path,'/',ids{i,1},'_V1_MR/gs_surf.mat']).gs_surf;
%             confounds(s,:)=load([path,'/',ids{i,1},'_V1_MR/confounds.mat']).confounds;
            temp1=load([path,'/',ids{i,1},'_V1_MR/censored_mask.mat']).censored_mask;
            temp2=load([path,'/',ids{i,1},'_V1_MR/Movement_RelativeRMS.mat']).Movement_RelativeRMS;
            Movement(s,1)=mean(temp2(temp1==1));
            break;
        end
    end
end

% partial correlation
coupling_age=zeros(atlas,4);
coupling_age_yeo=zeros(max(bna2yeo(:,1)),4);
for i=1:atlas
    mask = ~isnan(adjusted_r2(i,:)') & ~isinf(adjusted_r2(i,:)');
    [coupling_age(i,1),coupling_age(i,2)]=partialcorr(adjusted_r2(i,mask)',demography(mask,3),[demography(mask,1) tbv(mask,1) Movement(mask,1)]);
end
coupling_age(:,3)=mafdr(coupling_age(:,2),'BHFDR', true);
coupling_age(:,4)=coupling_age(:,2)*atlas;
save([savepath,'coupling_age.mat'],'coupling_age');
find(coupling_age(:,3)<0.05)
find(coupling_age(:,4)<0.05)
[temp1,temp2]=partialcorr(nanmean(adjusted_r2)',demography(:,3),[demography(mask,1) tbv(mask,1) Movement(mask,1)]);

coupling_age_yeo_anova=[0,0];
coupling_age_yeo_anova_pos=[0,0];
coupling_age_yeo_anova_neg=[0,0];
temp=0;
for i=1:max(bna2yeo(:,1))
    mask = ~isnan(adjusted_r2_yeo(i,:)') & ~isinf(adjusted_r2_yeo(i,:)');
    [coupling_age_yeo(i,1),coupling_age_yeo(i,2)]=partialcorr(adjusted_r2_yeo(i,mask)',demography(mask,3),[demography(mask,1) tbv(mask,1) Movement(mask,1)]);
    coupling_age_yeo_draw{i}=coupling_age(bna2yeo(:,1)==i,1);
    coupling_age_yeo_draw_pos{i}=coupling_age_yeo_draw{i}(coupling_age_yeo_draw{i}>0);
    coupling_age_yeo_draw_neg{i}=coupling_age_yeo_draw{i}(coupling_age_yeo_draw{i}<0);
    temp=temp+sum(coupling_age(bna2yeo(:,1)==i,1));
    coupling_age_yeo_anova_neg=[coupling_age_yeo_anova_neg;[coupling_age_yeo_draw_neg{i} ones(size(coupling_age_yeo_draw_neg{i}))*i]];
    coupling_age_yeo_anova_pos=[coupling_age_yeo_anova_pos;[coupling_age_yeo_draw_pos{i} ones(size(coupling_age_yeo_draw_pos{i}))*i]];
end
covariates=[demography(mask,1) tbv(mask,1) Movement(mask,1)];
save([savepath,'covariates.mat'],'covariates');

coupling_age_yeo(:,3)=mafdr(coupling_age_yeo(:,2),'BHFDR', true);
coupling_age_yeo(:,4)=coupling_age_yeo(:,2)*atlas;
save([savepath,'coupling_age_yeo.mat'],'coupling_age_yeo');

coupling_age_yeo_anova_pos(1,:)=[];
[~,tb1,stats1]=kruskalwallis(coupling_age_yeo_anova_pos(:,1),coupling_age_yeo_anova_pos(:,2),'off');
[post_hoc1,~,~,~]=multcompare(stats1);
coupling_age_yeo_anova_neg(1,:)=[];
[~,tb1,stats1]=kruskalwallis(coupling_age_yeo_anova_neg(:,1),coupling_age_yeo_anova_neg(:,2),'off');

coupling_age_yeo_anova(1,:)=[];
[~,tb1,stats1]=kruskalwallis(coupling_age_yeo_anova(:,1),coupling_age_yeo_anova(:,2),'off');
[post_hoc1,~,~,~]=multcompare(stats1);
[~,tb2,stats2]=anova1(coupling_age_yeo_anova(:,1),coupling_age_yeo_anova(:,2),'off');
[post_hoc2,~,~,~]=multcompare(stats2);

index=[6 2 3 7 4 1 5];
for i=1:7
    for j=1:7
        [~,p(i,j)]=ttest2(coupling_age_yeo_anova(coupling_age_yeo_anova(:,2)==index(i),1),coupling_age_yeo_anova(coupling_age_yeo_anova(:,2)==index(j),1));
    end
end


% GLM
slope_glm=load('D:\Work\hcp_development\results\coupling-multimodal\couplingwithage\glm\couplingslope.mat').slope;
slope=slope_glm(638,:);
slope_glm=slope_glm(1:atlas,:);
coupling_age_glm=zeros(size(slope_glm,1),5);
coupling_age_glm(:,1:3)=slope_glm(:,1:3);
coupling_age_glm(:,4)=mafdr(coupling_age_glm(:,3),'BHFDR', true);
coupling_age_glm(:,5)=coupling_age_glm(:,3)*atlas;
save([savepath,'coupling_age_glm.mat'],'coupling_age_glm');

coupling_age_glm_yeo_anova=[0,0];
temp=0;
for i=1:max(bna2yeo(:,1))
    coupling_age_glm_yeo_draw{i}=coupling_age_glm(bna2yeo(:,1)==i,1)*1000;
    coupling_age_glm_yeo_anova=[coupling_age_glm_yeo_anova;[coupling_age_glm_yeo_draw{i} ones(size(coupling_age_glm_yeo_draw{i}))*i]];
end

coupling_age_glm_yeo_anova(1,:)=[];

[~,tb1,stats1]=kruskalwallis(coupling_age_glm_yeo_anova(:,1),coupling_age_glm_yeo_anova(:,2),'off');
[post_hoc1,~,~,~]=multcompare(stats1);
[~,tb2,stats2]=anova1(coupling_age_glm_yeo_anova(:,1),coupling_age_glm_yeo_anova(:,2),'off');
[post_hoc2,~,~,~]=multcompare(stats2);
% for i=1:size(post_hoc2,1)
%     post_hoc2_p(post_hoc2(i,1),post_hoc2(i,2))=post_hoc2(i,6);
%     post_hoc2_p(post_hoc2(i,2),post_hoc2(i,1))=post_hoc2(i,6);
% end



coupling_age_glm_2=zeros(size(slope_glm,1),5);
coupling_age_glm_2(:,1)=slope_glm(:,10);%age2
coupling_age_glm_2(:,2)=-slope_glm(:,1)./slope_glm(:,10)*2;%peak
coupling_age_glm_2(:,3)=mafdr(slope_glm(:,12),'BHFDR', true);%age2 p
coupling_age_glm_2(:,4)=mafdr(slope_glm(:,5),'BHFDR', true);%anova p

save([savepath,'coupling_age_glm_2.mat'],'coupling_age_glm_2');
find(coupling_age_glm_2(:,3) & coupling_age_glm_2(:,4)<0.05)
coupling_age_glm_2(coupling_age_glm_2(:,3) & coupling_age_glm_2(:,4)<0.05,2)

adjusted_r2_mean_fitted=nanmean(adjusted_r2)'-slope(:,8).*demography(:,1)-slope(:,7).*tbv(:,1)-slope(:,6).*Movement(:,1);


