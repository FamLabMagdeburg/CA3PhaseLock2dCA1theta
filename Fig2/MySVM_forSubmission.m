% script MySVM_forSubmission
% calculate the performance using firing rates
% Also detect which neuron is effective by deleting neuron one by one.
% Matlab script for the paper "Phase-locking of hippocampal CA3 neurons to
% distal CA1 theta oscillations selectively predicts memory performance "
% Fig.2A-C
% start
% Go to the folder where the parameter and spiking data are saved
load animal_params.mat
% % all data
file_tet=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; %the list denotes to the numbering of the tetrodes where the spikes were recorded
data=f_extract_spikes_Neuralynx('animal_tet',file_tet,000, 2000,ts_stimon_exp,0); % extracting the 2000ms firing data for SVM analysis  

% extract features
ii=0;
neuron_id=[];
for i=1:length(data)
    for j=1:length(data(i).ss)
        ii=ii+1;
        neuron_id(1,ii)=file_tet(i);
        neuron_id(2,ii)=j;
        for k=1:length(data(i).ss(j).xx)
            meanfr(k,ii)=length(data(i).ss(j).xx(k).times);
        end
    end
end

%
trial=1:length(label_oldnew);
label=label_oldnew; %labels of different trial types, old=repeated trials; new=non-repeated trials
meanfrM=meanfr;


for i=1:length(trial)
label_train=label';
label_train(i)=[];
meanfr_train=meanfrM;
meanfr_train(i,:)=[];
model = svmtrain(label_train, meanfr_train, '-t 1');
[predicted_label] = svmpredict(label(i), squeeze(meanfrM(i,:)), model);
test2(i)=predicted_label;
end
perf_all=length(find(test2==label))/length(trial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the effect of deleting any neuron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check_meanfr=meanfr;
perf_effected=[];
for j=1:size(check_meanfr,2)
    meanfrM=check_meanfr;
    meanfrM(:,j)=[];
    meanfrM=squeeze(meanfrM);
    for i=1:length(trial)
    label_train=label';
    label_train(i)=[];
    label_train=squeeze(label_train);
    meanfr_train=meanfrM;
    meanfr_train(i,:)=[];
    meanfr_train=squeeze(meanfr_train);
    model = svmtrain(label_train, meanfr_train, '-s 0 -t 1');
    [predicted_label , accuracy, decision_values] = svmpredict(label(i), squeeze(meanfrM(i,:)), model);
        test2(i)=predicted_label;
end
perf=length(find(test2==label))/length(trial);
perf_effected(j)=perf-perf_all;
end
% Are the effective neurons enough for the high performance
xx=find(perf_effected<0);
meanfrM=meanfr(:,xx);
for i=1:length(trial)
label_train=label';
label_train(i)=[];
meanfr_train=meanfrM;
meanfr_train(i,:)=[];
% libsvm
model = svmtrain(label_train, meanfr_train, '-t 1');
[predicted_label] = svmpredict(label(i), squeeze(meanfrM(i,:)), model);

% % Matlab
% model = fitcsvm( meanfr_train,label_train,'KernelFunction','rbf','Standardize',true)
% [predicted_label] = predict( model,squeeze(meanfrM(i,:)));

test2(i)=predicted_label;
end
perf_effec=length(find(test2==label))/length(trial)
% How is the performance of the neurons without the effective ones
meanfrM=meanfr;
meanfrM(:,xx)=[];
meanfrM=squeeze(meanfrM);
for i=1:length(trial)
label_train=label';
label_train(i)=[];
meanfr_train=meanfrM;
meanfr_train(i,:)=[];
% libsvm
model = svmtrain(label_train, meanfr_train, '-t 1');
[predicted_label] = svmpredict(label(i), squeeze(meanfrM(i,:)), model);

% % Matlab
% model = fitcsvm( meanfr_train,label_train,'KernelFunction','rbf','Standardize',true)
% [predicted_label] = predict( model,squeeze(meanfrM(i,:)));

test2(i)=predicted_label;
end
perf_others=length(find(test2==label))/length(trial)

