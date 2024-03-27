% A script for analyzing the moving speed of one example animal.
% Including extracting data from position files (fter DeepLabCut), convert to speed
% (cm/sec), extracting theta amplitudes cycle-by-cycle, then compare their
% rlationships.
% For the paper titled "Phase-locking of hippocampal CA3 neurons to distal
% CA1 theta oscillations selectively predicts memory performance. "
% Author Shihpi Ku
% date 13.10.2022

%% example animal LE84
clear all
close all
animal='LE84';
x_ppc=405.01/35;    %pixel/cm
y_ppc=207/20;   %pixel/cm

%%
%cd \Data\DEEPLABCUT\videos
% study phase
M84_study=csvread('\\DEEPLABCUT\Video_cuts_timestamp\vlc-record-2022-10-04-16h59m37s-compressed_LE84_20190712_exp_studyDLC_resnet50_LE_4ratsSep22shuffle1_200000.csv',3,0);
ears84_study(:,1)=(M84_study(:,8)+M84_study(:,11))./2;
ears84_study(:,2)=(M84_study(:,9)+M84_study(:,12))./2;
ears84_study(:,3)=M84_study(:,1).*0.01+0.01;
[t84_study,x84_study,y84_study,vx84_study,vy84_study,ax84_study,ay84_study] = KalmanVel(ears84_study(:,1),ears84_study(:,2),ears84_study(:,3),2);
x84_study=x84_study/x_ppc;
y84_study=y84_study/y_ppc;

% [h,p]=ttest2(vy84_bsl,vy84_test)
% [h,p]=ttest2(vx84_bsl,vx84_test)
LE84_study_t=round(xlsread('\\DEEPLABCUT\Proximodistal_Video_TimeStamps.xlsx','LE84_20190712','E41:E50').*100)-129279;    %129279: the starting time stamp of study phase
% Get all the study trials, 2000ms per trial, in total 200 frames (10ms/frame)
for i=1:10
    M84_study_trial((i-1)*200+1:i*200,1,1)=vx84_study(LE84_study_t(i)-200:LE84_study_t(i)-1);   
    M84_study_trial((i-1)*200+1:i*200,2,1)=vy84_study(LE84_study_t(i):LE84_study_t(i)+199);
    M84_study_trial((i-1)*200+1:i*200,1,2)=vx84_study(LE84_study_t(i)+200:LE84_study_t(i)+399);
    M84_study_trial((i-1)*200+1:i*200,2,2)=vy84_study(LE84_study_t(i)-200:LE84_study_t(i)-1);
    M84_study_trial((i-1)*200+1:i*200,1,3)=vx84_study(LE84_study_t(i):LE84_study_t(i)+199);
    M84_study_trial((i-1)*200+1:i*200,2,3)=vy84_study(LE84_study_t(i)+200:LE84_study_t(i)+399);
    for j=1:200
        temp=sqrt((vx84_study(round(LE84_study_t(i))+j-1)/x_ppc)^2+(vy84_study(round(LE84_study_t(i))+j-1)/y_ppc)^2);
        dist84_study(i,j)=temp;
    end
end

%% baseline
M84_bsl = csvread('\\DEEPLABCUT\Video_cuts_timestamp\vlc-record-2022-10-04-16h53m22s-compressed_LE84_20190712_exp_bslDLC_resnet50_LE_Proximodistal_bslSep23shuffle1_200000.csv',3,0);
ears84_bsl(:,1)=(M84_bsl(:,8)+M84_bsl(:,11))./2;
ears84_bsl(:,2)=(M84_bsl(:,9)+M84_bsl(:,12))./2;
ears84_bsl(:,3)=M84_bsl(:,1).*0.01+0.01;
[t84_bsl,x84_bsl,y84_bsl,vx84_bsl,vy84_bsl,ax84_bsl,ay84_bsl] = KalmanVel(ears84_bsl(:,1),ears84_bsl(:,2),ears84_bsl(:,3),2);
x84_bsl=x84_bsl/x_ppc;
y84_bsl=y84_bsl/y_ppc;
LE84_bsl_t=LE84_study_t;
dist84_bsl=[];
for i=1:10
    for j=1:200
        temp=sqrt((vx84_bsl(round(LE84_bsl_t(i))+j-1)/x_ppc)^2+(vy84_bsl(round(LE84_bsl_t(i))+j-1)/y_ppc)^2);
        dist84_bsl(i,j)=temp;
    end
end
%% test phase
M84_test=csvread('\\DEEPLABCUT\Video_cuts_timestamp\vlc-record-2022-10-04-17h05m31s-compressed_LE84_20190712_exp_testDLC_resnet50_LE_4ratsSep22shuffle1_200000.csv',3,0);
ears84_test(:,1)=(M84_test(:,8)+M84_test(:,11))./2;
ears84_test(:,2)=(M84_test(:,9)+M84_test(:,12))./2;
ears84_test(:,3)=M84_test(:,1).*0.01+0.01;
[t84_test,x84_test,y84_test,vx84_test,vy84_test,ax84_test,ay84_test] = KalmanVel(ears84_test(:,1),ears84_test(:,2),ears84_test(:,3),2);
x84_test=x84_test/x_ppc;
y84_test=y84_test/y_ppc;

LE84_test_t=round(xlsread('\\DEEPLABCUT\Proximodistal_Video_TimeStamps.xlsx','LE84_20190712','L41:L60').*100)-274333;
for i=1:20
    M84_test_trial((i-1)*200+1:i*200,1,1)=vx84_test(LE84_test_t(i)-200:LE84_test_t(i)-1);
    M84_test_trial((i-1)*200+1:i*200,2,1)=vy84_test(LE84_test_t(i):LE84_test_t(i)+199);
    M84_test_trial((i-1)*200+1:i*200,1,2)=vx84_test(LE84_test_t(i)+200:LE84_test_t(i)+399);
    M84_test_trial((i-1)*200+1:i*200,2,2)=vy84_test(LE84_test_t(i)-200:LE84_test_t(i)-1);
    M84_test_trial((i-1)*200+1:i*200,1,3)=vx84_test(LE84_test_t(i):LE84_test_t(i)+199);
    M84_test_trial((i-1)*200+1:i*200,2,3)=vy84_test(LE84_test_t(i)+200:LE84_test_t(i)+399);
    for j=1:200
        temp=sqrt((vx84_test(round(LE84_test_t(i))+j-1)/x_ppc)^2+(vy84_test(round(LE84_test_t(i))+j-1)/y_ppc)^2);
        dist84_test(i,j)=temp;
    end
end

load \\Group_analysis\LE84_20190917_params label_oldnew
old_trial=find(label_oldnew==0);
new_trial=find(label_oldnew==1);
dist84_old=dist84_test(old_trial,:);
dist84_new=dist84_test(new_trial,:);
%%
[h,p]=ttest2(vy84_study,vy84_test)
[h,p]=ttest2(vx84_study,vx84_test)
[h,p]=ttest2(ax84_study,ax84_test)
[h,p]=ttest2(ay84_study,ay84_test)

%% linear fit and plot
% figure
% plot(study_M(1,:),study_M(2,:),'.b')
% hold on
% yhat=predict(mdStudy,[0;45000]);
% plot([0;45000],yhat,'-k')

%% read theta cycles

%%
% study phase
theta_M_study.data=[];
theta_M_study.dataname{1}='trial_no';
theta_M_study.dataname{2}='theta_amp';
theta_M_study.dataname{3}='cycle_time';
theta_M_study.dataname{4}='speed';

ii=0;
for i=1:10
    data=readtable('\\data\theta_test_LE84.xlsx','sheet',['LE84_2_' num2str(i)]);
    theta_M_study.data(ii+1:ii+length(data.volt_amp),1)=ones(length(data.volt_amp),1).*i;
    theta_M_study.data(ii+1:ii+length(data.volt_amp),2)=data.volt_amp;
    theta_M_study.data(ii+1:ii+length(data.volt_amp),3)=round(((data.sample_last_trough+data.sample_next_trough)/2)/40);
    theta_M_study.data(ii+1:ii+length(data.volt_amp),4)=0;
    for j=1:+length(data.volt_amp)
        theta_M_study.data(ii+j,4)=mean(dist84_study(i,theta_M_study.data(ii+j,3)-4:theta_M_study.data(ii+j,3)+3));
    end
    ii=ii+length(data.volt_amp);
end
%%
% baseline
theta_M_bsl.data=[];
theta_M_bsl.dataname{1}='trial_no';
theta_M_bsl.dataname{2}='theta_amp';
theta_M_bsl.dataname{3}='cycle_time';
theta_M_bsl.dataname{4}='speed';

ii=0;
for i=1:10
    data=readtable('\\data\theta_test_LE84.xlsx','sheet',['LE84_1_' num2str(i)]);
    theta_M_bsl.data(ii+1:ii+length(data.volt_amp),1)=ones(length(data.volt_amp),1).*i;
    theta_M_bsl.data(ii+1:ii+length(data.volt_amp),2)=data.volt_amp;
    theta_M_bsl.data(ii+1:ii+length(data.volt_amp),3)=round(((data.sample_last_trough+data.sample_next_trough)/2)/40);
    theta_M_bsl.data(ii+1:ii+length(data.volt_amp),4)=0;
    for j=1:+length(data.volt_amp)
        theta_M_bsl.data(ii+j,4)=mean(dist84_bsl(i,theta_M_bsl.data(ii+j,3)-4:theta_M_bsl.data(ii+j,3)+3));
    end
    ii=ii+length(data.volt_amp);
end
% test phase
theta_M_old.data=[];
theta_M_old.dataname{1}='trial_no';
theta_M_old.dataname{2}='theta_amp';
theta_M_old.dataname{3}='cycle_time';
theta_M_old.dataname{4}='speed';

ii=0;
% test old
for i=1:10
    data=readtable('\\data\theta_test_LE84.xlsx','sheet',['LE84_3_' num2str(i)]);
    theta_M_old.data(ii+1:ii+length(data.volt_amp),1)=ones(length(data.volt_amp),1).*i;
    theta_M_old.data(ii+1:ii+length(data.volt_amp),2)=data.volt_amp;
    theta_M_old.data(ii+1:ii+length(data.volt_amp),3)=round(((data.sample_last_trough+data.sample_next_trough)/2)/40);
    theta_M_old.data(ii+1:ii+length(data.volt_amp),4)=0;
    for j=1:+length(data.volt_amp)
        theta_M_old.data(ii+j,4)=mean(dist84_old(i,theta_M_old.data(ii+j,3)-4:theta_M_old.data(ii+j,3)+3));
    end
    ii=ii+length(data.volt_amp);
end
% test new
theta_M_new.data=[];
theta_M_new.dataname{1}='trial_no';
theta_M_new.dataname{2}='theta_amp';
theta_M_new.dataname{3}='cycle_time';
theta_M_new.dataname{4}='speed';

ii=0;

for i=1:10
    data=readtable('\\data\theta_test_LE84.xlsx','sheet',['LE84_4_' num2str(i)]);
    theta_M_new.data(ii+1:ii+length(data.volt_amp),1)=ones(length(data.volt_amp),1).*i;
    theta_M_new.data(ii+1:ii+length(data.volt_amp),2)=data.volt_amp;
    theta_M_new.data(ii+1:ii+length(data.volt_amp),3)=round(((data.sample_last_trough+data.sample_next_trough)/2)/40);
    theta_M_new.data(ii+1:ii+length(data.volt_amp),4)=0;
    for j=1:+length(data.volt_amp)
        theta_M_new.data(ii+j,4)=mean(dist84_new(i,theta_M_new.data(ii+j,3)-4:theta_M_new.data(ii+j,3)+3));
    end
    ii=ii+length(data.volt_amp);
end
theta_M_test=[theta_M_old.data;theta_M_new.data];

%%
figure
plot(theta_M_old.data(:,4),theta_M_old.data(:,2),'.r')
hold on
plot(theta_M_new.data(:,4),theta_M_new.data(:,2),'.b')

plot(theta_M_study.data(:,4),theta_M_study.data(:,2),'.k')

mdl_s = fitlm(theta_M_study.data(find(theta_M_study.data(:,4)~=0),4),theta_M_study.data(find(theta_M_study.data(:,4)~=0),2))
mdl_o = fitlm(theta_M_old.data(find(theta_M_old.data(:,4)~=0),4),theta_M_old.data(find(theta_M_old.data(:,4)~=0),2))
mdl_n = fitlm(theta_M_new.data(find(theta_M_new.data(:,4)~=0),4),theta_M_new.data(find(theta_M_new.data(:,4)~=0),2))
mdl_t = fitlm(theta_M_test(:,4),theta_M_test(:,2))
mdl_b = fitlm(theta_M_bsl.data(find(theta_M_bsl.data(:,4)~=0),4),theta_M_bsl.data(find(theta_M_bsl.data(:,4)~=0),2))
%% plot study vs. test
figure
h=plot(mdl_s)
% change plot color according to https://de.mathworks.com/matlabcentral/answers/538543-how-to-change-the-color-of-different-parts-of-regression-plots-using-fitlm-function
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
% The confidence bounds have 2 handles but only one of 
% the handles contains the legend string.  The first
% line below finds that object and then searches for 
% other objects in the plot that have the same linestyle
% and color. 
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'k';
set(cbHandles, 'Color', 'k', 'LineWidth', 1)
fitHandle.Color = [0 0 0]; %black, or fitHandle.Color = [1 .6445 0]; %orange 
hold on

%plot(mdl_o,'Color','r')
ht=plot(mdl_t,'Color','r')
xlabel('speed (cm/s)')
ylabel('theta amplitude')
title(animal)
legend([h(1) ht(1)],{'study','test'})
%axis([0 400 2000 16000])
%% plot baseline, study vs. test
figure
h=plot(mdl_s)
% change plot color according to https://de.mathworks.com/matlabcentral/answers/538543-how-to-change-the-color-of-different-parts-of-regression-plots-using-fitlm-function
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
% The confidence bounds have 2 handles but only one of 
% the handles contains the legend string.  The first
% line below finds that object and then searches for 
% other objects in the plot that have the same linestyle
% and color. 
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'k';
set(cbHandles, 'Color', 'k', 'LineWidth', 1)
fitHandle.Color = [0 0 0]; %black, or fitHandle.Color = [1 .6445 0]; %orange 
hold on
hb=plot(mdl_b)
% change plot color according to https://de.mathworks.com/matlabcentral/answers/538543-how-to-change-the-color-of-different-parts-of-regression-plots-using-fitlm-function
% Get handles to plot components
dataHandle = findobj(hb,'DisplayName','Data');
fitHandle = findobj(hb,'DisplayName','Fit');
% The confidence bounds have 2 handles but only one of 
% the handles contains the legend string.  The first
% line below finds that object and then searches for 
% other objects in the plot that have the same linestyle
% and color. 
cbHandles = findobj(hb,'DisplayName','Confidence bounds');
cbHandles = findobj(hb,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'b';
set(cbHandles, 'Color', 'b', 'LineWidth', 1)
fitHandle.Color = [0 0 1]; %black, or fitHandle.Color = [1 .6445 0]; %orange 

%plot(mdl_o,'Color','r')
ht=plot(mdl_t,'Color','r')
xlabel('speed (cm/s)')
ylabel('theta amplitude')
title(animal)
legend([hb(1) h(1) ht(1)],{'baseline','study','test'})
%% plot test repeated vs. non-repeated
figure
h=plot(mdl_n)
% change plot color according to https://de.mathworks.com/matlabcentral/answers/538543-how-to-change-the-color-of-different-parts-of-regression-plots-using-fitlm-function
% Get handles to plot components
dataHandle = findobj(h,'DisplayName','Data');
fitHandle = findobj(h,'DisplayName','Fit');
% The confidence bounds have 2 handles but only one of 
% the handles contains the legend string.  The first
% line below finds that object and then searches for 
% other objects in the plot that have the same linestyle
% and color. 
cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
dataHandle.Color = 'g';
set(cbHandles, 'Color', 'g', 'LineWidth', 1)
fitHandle.Color = [0 1 0]; %black, or fitHandle.Color = [1 .6445 0]; %orange 
hold on
plot(mdl_o,'Color','r')
xlabel('speed(cm/s)')
ylabel('theta amplitude')
title([animal])
legend off
%axis([0 400 2000 16000])

%% plot the histogram of speed and theta amp
% speed
% smooth the speed, averaged every 80ms (8 frames), no overlap
for i=1:25
dist84s_bsl(:,i)=mean(dist84_bsl(:,(i-1)*8+1:i*8),2);
end
for i=1:25
dist84s_study(:,i)=mean(dist84_study(:,(i-1)*8+1:i*8),2);
end
for i=1:25
dist84s_test(:,i)=mean(dist84_test(:,(i-1)*8+1:i*8),2);
end
dist84s_old=dist84s_test(old_trial,:);
dist84s_new=dist84s_test(new_trial,:);
figure
boxplot([dist84s_bsl(:) dist84s_study(:) dist84s_old(:) dist84s_new(:)])
xticklabels({'baseline','study','repeated','non-repeated'})
ylabel('speed (cm/sec)')
title([animal ' speed'])
%%
figure
hb=histfit(dist84s_bsl(:),100)
hold on
hs=histfit(dist84s_study(:),100)
ht=histfit(dist84s_test(:),100)
alpha(.5)
hs(1).FaceColor=[0.2 0.2 .2];
hs(2).Color=[0 0 0];
hb(2).Color=[0 0 1];
hb(1).EdgeColor='none';
legend([hb(1) hs(1) ht(1)],{'baseline','study','test'})
title([animal ' speed density'])
xlabel('speed(cm/s)')
ylabel('count')

%% theta amplitude
figure
ht_b=histfit(theta_M_bsl.data(:,2),100)
hold on
ht_s=histfit(theta_M_study.data(:,2),100)
ht_t=histfit(theta_M_test(:,2),100)
alpha(.5)
ht_s(1).FaceColor=[0.2 0.2 .2];
ht_s(2).Color=[0 0 0];
ht_b(2).Color=[0 0 1];
legend([ht_b(1) ht_s(1) ht_t(1)],{'baseline','study','test'})
title([animal ' theta amplitude density'])
xlabel('theta amplitude')
ylabel('count')