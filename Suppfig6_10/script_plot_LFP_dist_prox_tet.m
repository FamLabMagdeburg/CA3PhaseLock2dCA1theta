% For the paper titled "Phase-locking of hippocampal CA3 neurons to distal CA1 theta oscillations selectively predicts memory performance" 
% Plot example LFP traces for reviewer's additional comment to show proxCA1 theta is strong and detectable. Overlay distCA1, seg2 and proximal CA1 traces. Need to manually insert the tetrode number for distal seg2, proximal tetrode for each animal. Did only LE82, still need to insert for other animals.   
% Now as Supplementary Fig. 10

clear all
close all

cd('R:\Data\forNico\LE84\LE84_20190712')
load LE84_20190917_expRAW_params.mat
animal='LE84_20190917_tet';
animalname='LE84';
bsl=2000;    % for raster plot
post=2000;  % for raster plot
bsl_lfp=000;     % for LFP
post_lfp=2000;     % for LFP
LFP_filename_pre='C:\Users\sku\LIN\data\Neuralynx\LE84_2019-07-12_13-15-11';
LFP_filename_post='_0001';
region_id=3;
tet_dist=2; 
tet_prox=12; 
trial=2;
seg2=7;

%%
ts=ts_stimon_exp;
ts_bsl=ts_stimon_study(1)-5000:-10000:ts_stimon_study(1)-200000;
trial=1:length(label_oldnew);
label=label_oldnew;
label_old=find(label==0);
label_new=find(label==1);
%params.tapers=[3 5];
params.tapers=[5 9]; 
params.err=[1 0.05]; % population error bars
params.Fs=4000;
params.trialave=1;
params.fpass=[1 135];
movingwin =[0.5 0.05];    % for theta

lfp_length=round((params.Fs/1000)*(bsl_lfp+post_lfp));
filter_band=300;
Nlx_samplerate=32000;
downsample_rate=8;
[b,a] = butter(5,filter_band./(Nlx_samplerate/2),'low');

%% extract dist tet LFP
ii=0;
data_dist={};
z=tet_dist;
ii=ii+1;
%filename=[animal num2str(tet) '.mat']
% extract LFP
[Samples_dist, Header] =f_extract_LFP_Nlx(LFP_filename_pre,LFP_filename_post, tet_dist);

lfp_dist_raw=Samples_dist(:);
disp(['low pass lfp filtering under ' num2str(filter_band) 'Hz....'])
lfp_dist_f=filtfilt(b,a,double(Samples_dist(:)));
lfp_dist=downsample(lfp_dist_f,downsample_rate,4);   %downsample 8 times so the Fs down from 32K to 4K, take the 4th element.
datalfp_dist=[];
datalfp_dist_study=[];
datalfp_dist_bsl=[];
% cut into trials
for i=1:length(ts)
    temp=double(extractdatac(lfp_dist,params.Fs,[ts(i)-bsl_lfp ts(i)+post_lfp]./1000));
    temp_bsl=double(extractdatac(lfp_dist,params.Fs,[ts_bsl(i)-bsl_lfp ts_bsl(i)+post_lfp]./1000));
    datalfp_dist(:,i)=temp(1:lfp_length)';
    datalfp_dist_bsl(:,i)=temp_bsl(1:lfp_length)';
end
%
for i=1:length(ts_stimon_study)
    temp_study=double(extractdatac(lfp_dist,params.Fs,[ts_stimon_study(i)-bsl_lfp ts_stimon_study(i)+post_lfp]./1000));
    datalfp_dist_study(:,i)=temp_bsl(1:lfp_length)';
end


%% look at the power spectrum during a trial or study and test phases.
     [S_dist_bsl,f,Serr_dist_bsl] = mtspectrumc( datalfp_dist_bsl,  params );
    [S_dist_study,f,Serr_dist_study] = mtspectrumc( datalfp_dist_study, params );
    [S_dist_new,f,Serr_dist_new] = mtspectrumc( datalfp_dist(:,label_new),  params );
    [S_dist_old,f,Serr_dist_old] = mtspectrumc( datalfp_dist(:,label_old), params );
%%   
figure
subplot(1,2,1)
plot(f,(S_dist_new),'g')
hold on
plot(f,(S_dist_old),'r')
plot(f,(S_dist_bsl),'k')
plot(f,(S_dist_study),'b')
title('dist')
subplot(1,2,2)
lower_band=find(f<=15);
plot(f(lower_band),S_dist_new(lower_band),'g')
hold on
plot(f(lower_band),S_dist_old(lower_band),'r')
plot(f(lower_band),S_dist_bsl(lower_band),'k')
plot(f(lower_band),S_dist_study(lower_band),'b')
title('zoom in lower band')
%%
lfp_dist_long=double(extractdatac(lfp_dist,params.Fs,[ts(1)-bsl_lfp ts(end)+post_lfp]./1000));
[S_dist,f_dist,Serr_dist] = mtspectrumc(lfp_dist_long,  params );
figure
subplot(1,2,1)
plot(f_dist,S_dist)
title('dist')
subplot(1,2,2)
lower_band=find(f_dist<=15);
plot(f_dist(lower_band),S_dist(lower_band))
delta_band=find(f_dist>=2 & f_dist<=4)
theta_band=find(f_dist>=6 & f_dist<=12)
ratio_dist=mean(S_dist(theta_band))/mean(S_dist(delta_band))
%%
lfp_dist_bsl_long=double(extractdatac(lfp_dist,params.Fs,[min(ts_bsl) max(ts_bsl)]./1000));
[S_dist_bsl,f_dist_bsl,Serr_dist_bsl] = mtspectrumc(lfp_dist_bsl_long,  params );
figure
subplot(1,2,1)
plot(f_dist_bsl,S_dist_bsl)
subplot(1,2,2)
lower_band=find(f_dist_bsl<=15);
plot(f_dist_bsl(lower_band),S_dist_bsl(lower_band))
delta_band=find(f_dist_bsl>=2 & f_dist_bsl<=4);
theta_band=find(f_dist_bsl>=6 & f_dist_bsl<=12);
ratio_dist_bsl=mean(S_dist_bsl(theta_band))/mean(S_dist_bsl(delta_band))
title('dist baseline')

%% extract prox tet LFP
ii=0;
data_prox={};
z=tet_prox;
ii=ii+1;
% extract LFP
[Samples_prox, Header] =f_extract_LFP_Nlx(LFP_filename_pre,LFP_filename_post, tet_prox);

lfp_prox_raw=Samples_prox(:);
disp(['low pass lfp filtering under ' num2str(filter_band) 'Hz....'])
lfp_prox_f=filtfilt(b,a,double(Samples_prox(:)));
lfp_prox=downsample(lfp_prox_f,downsample_rate,4);   %downsample 8 times so the Fs down from 32K to 4K, take the 4th element.
datalfp_prox=[];
datalfp_prox_study=[];
datalfp_prox_bsl=[];
% cut into trials
for i=1:length(ts)
    temp=double(extractdatac(lfp_prox,params.Fs,[ts(i)-bsl_lfp ts(i)+post_lfp]./1000));
    temp_bsl=double(extractdatac(lfp_prox,params.Fs,[ts_bsl(i)-bsl_lfp ts_bsl(i)+post_lfp]./1000));
    
    %     if corerr(i)>0
    %         ii=ii+1;
    datalfp_prox(:,i)=temp(1:lfp_length)';
    datalfp_prox_bsl(:,i)=temp_bsl(1:lfp_length)';
    %     end
end
%
for i=1:length(ts_stimon_study)
    temp_study=double(extractdatac(lfp_prox,params.Fs,[ts_stimon_study(i)-bsl_lfp ts_stimon_study(i)+post_lfp]./1000));
    datalfp_prox_study(:,i)=temp_bsl(1:lfp_length)';
end


     [S_prox_bsl,f,Serr_prox_bsl] = mtspectrumc( datalfp_prox_bsl,  params );
    [S_prox_study,f,Serr_prox_study] = mtspectrumc( datalfp_prox_study, params );
    [S_prox_new,f,Serr_prox_new] = mtspectrumc( datalfp_prox(:,label_new),  params );
    [S_prox_old,f,Serr_prox_old] = mtspectrumc( datalfp_prox(:,label_old), params );
%%
%% extract seg2 tet LFP

ii=0;
data_seg2={};
z=seg2;
ii=ii+1;
% extract LFP
[Samples_seg2, Header] =f_extract_LFP_Nlx(LFP_filename_pre,LFP_filename_post, seg2);

lfp_seg2_raw=Samples_seg2(:);
disp(['low pass lfp filtering under ' num2str(filter_band) 'Hz....'])
lfp_seg2_f=filtfilt(b,a,double(Samples_seg2(:)));
lfp_seg2=downsample(lfp_seg2_f,downsample_rate,4);   %downsample 8 times so the Fs down from 32K to 4K, take the 4th element.
datalfp_seg2=[];
datalfp_seg2_study=[];
datalfp_seg2_bsl=[];
% cut into trials
for i=1:length(ts)
    temp=double(extractdatac(lfp_seg2,params.Fs,[ts(i)-bsl_lfp ts(i)+post_lfp]./1000));
    temp_bsl=double(extractdatac(lfp_seg2,params.Fs,[ts_bsl(i)-bsl_lfp ts_bsl(i)+post_lfp]./1000));
    datalfp_seg2(:,i)=temp(1:lfp_length)';
    datalfp_seg2_bsl(:,i)=temp_bsl(1:lfp_length)';
end
%
for i=1:length(ts_stimon_study)
    temp_study=double(extractdatac(lfp_seg2,params.Fs,[ts_stimon_study(i)-bsl_lfp ts_stimon_study(i)+post_lfp]./1000));
    datalfp_seg2_study(:,i)=temp_bsl(1:lfp_length)';
end
%%
figure
subplot(1,2,1)
plot(f,(S_prox_new),'g')
hold on
plot(f,(S_prox_old),'r')
plot(f,(S_prox_bsl),'k')
plot(f,(S_prox_study),'b')
title('prox')

subplot(1,2,2)
lower_band=find(f<=15);
plot(f(lower_band),S_prox_new(lower_band),'g')
hold on
plot(f(lower_band),S_prox_old(lower_band),'r')
plot(f(lower_band),S_prox_bsl(lower_band),'k')
plot(f(lower_band),S_prox_study(lower_band),'b')
title('zoom in lower band')
%%

lfp_prox_long=double(extractdatac(lfp_prox,params.Fs,[ts(1)-bsl_lfp ts(end)+post_lfp]./1000));
[S_prox,f_prox,Serr_prox] = mtspectrumc(lfp_prox_long,  params );
figure
subplot(1,2,1)
plot(f_prox,S_prox)
title('prox')

subplot(1,2,2)
lower_band=find(f_prox<=15);
plot(f_prox(lower_band),S_prox(lower_band))
delta_band=find(f>=2 & f<=4)
theta_band=find(f>=6 & f<=12)
ratio_prox=mean(S_prox(theta_band))/mean(S_prox(delta_band))


tstamp=linspace(0,size(datalfp_dist,1)/params.Fs,size(datalfp_dist,1));
Fs=params.Fs;
%%
trial=label_old(1);
figure
clf
hold on
plot(tstamp,bandpass(datalfp_prox(:,trial),[6 12],params.Fs),'k','LineWidth',0.5);
plot(tstamp,bandpass(datalfp_dist(:,trial),[6 12],params.Fs),'r','LineWidth',0.5);
legend('seg4(prox)','seg1(dist)')
title([animalname ' trial ' num2str(trial)])
xlabel('time(sec)')

%
figure
clf
tstamp=linspace(0,size(datalfp_dist,1)/params.Fs,size(datalfp_dist,1));
hold on
plot(tstamp,highpass(datalfp_prox(:,trial),4,params.Fs),'Color',[0.5 0.5 0.5],'LineWidth',0.5);
plot(tstamp,highpass(datalfp_dist(:,trial),4,params.Fs),'r','LineWidth',0.5);
legend('seg4(prox)','seg1(dist)')
title([animalname ' trial ' num2str(trial)])
xlabel('time(sec)')
figure
clf
hold on
plot(tstamp,bandpass(datalfp_prox(:,trial),[6 12],params.Fs),'k','LineWidth',0.5);
plot(tstamp,bandpass(datalfp_seg2(:,trial),[6 12],params.Fs),'b','LineWidth',0.5);
plot(tstamp,bandpass(datalfp_dist(:,trial),[6 12],params.Fs),'r','LineWidth',0.5);

legend('seg4(prox)','seg2','seg1(dist)')
title([animalname ' trial ' num2str(trial)])
xlabel('time(sec)')

%
figure
clf
tstamp=linspace(0,size(datalfp_dist,1)/params.Fs,size(datalfp_dist,1));
hold on
plot(tstamp,highpass(datalfp_prox(:,trial),4,params.Fs),'Color',[0.5 0.5 0.5],'LineWidth',0.5);

plot(tstamp,highpass(datalfp_seg2(:,trial),4,params.Fs),'b','LineWidth',0.5);
plot(tstamp,highpass(datalfp_dist(:,trial),4,params.Fs),'r','LineWidth',0.5);
legend('seg4(prox)','seg2','seg1(dist)')
title([animalname ' trial ' num2str(trial)])
xlabel('time(sec)')


