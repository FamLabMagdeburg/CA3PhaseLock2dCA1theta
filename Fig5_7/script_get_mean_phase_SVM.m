% get the theta phase of each individual animal and perform SVM 

clear all
close all
load('LE82_coherence_theta_all_phase.mat')  
animal_ID=1;
params='-s 0 -t 1 -c 10';


% 
phase_CA3_CA3=[];
phase_CA3_CA1=[];
phase_CA1_CA1=[];
phase_CA3_CA3_old_all={};
count_CA3_CA1_old=zeros(1,10);
phase_CA1_CA1_new_all={};
count_CA1_CA1_new=zeros(1,10);
thr=0.05;
load('region.mat');
CA3=find(region(animal_ID).area(:,1)==3)
CA1=find(region(animal_ID).area(:,1)==1)
% use dCA1 as reference theta
if ~isempty(find(region(animal_ID).area(:,1)==1 & region(animal_ID).area(:,2)==1))
    CA1_ref=find(region(animal_ID).area(:,1)==1 & region(animal_ID).area(:,2)==1);
else
    CA1_ref=find(region(animal_ID).area(:,1)==1 & region(animal_ID).area(:,2)==2);
end

% use pCA1 as reference theta
% if ~isempty(find(region(animal_ID).area(:,1)==1 & region(animal_ID).area(:,2)==5))
%     CA1_ref=find(region(animal_ID).area(:,1)==1 & region(animal_ID).area(:,2)==5);
% else
%     CA1_ref=find(region(animal_ID).area(:,1)==1 & region(animal_ID).area(:,2)==4);
% end

ii=0;
jj=0;
kk=0;

for i=1:length(data)
    % get CA1_CA1 phase
    if find(CA1==data{i}.tet_spk)
        for m=1:length(data{i}.coherence{data{i}.tet_spk}.neuron)
            for tt=1:10
                if ~isnan(mean(data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_new(:,tt)))
                    count_CA1_CA1_new(1,tt)=count_CA1_CA1_new(1,tt)+1;
                    phase_CA1_CA1_new_all{tt}(count_CA1_CA1_new(1,tt))=circ_mean(data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_new(:,tt));
                end
            end
            
            if isempty(find(isnan(data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_old))) && isempty(find(isnan(data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_new)))
                ii=ii+1;
                phase_CA1_CA1(ii,:)=[circ_mean((data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_old)) circ_mean((data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_new))];
            end
        end
    end
    % get CA3-CA3 phase
    if find(CA3==data{i}.tet_spk)
        for m=1:length(data{i}.coherence{data{i}.tet_spk}.neuron)
            if isempty(find(isnan(data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_old))) && isempty(find(isnan(data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_new)))
                jj=jj+1;
                phase_CA3_CA3(jj,:)=[circ_mean((data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_old)) circ_mean((data{i}.coherence{data{i}.tet_spk}.neuron(m).phi_new))];
            end
        end
        % get CA3-CA1 phase
                for m=1:length(data{i}.coherence{CA1_ref(1)}.neuron)
                    for tt=1:10
                        if ~isnan(mean(data{i}.coherence{CA1_ref(1)}.neuron(m).phi_old(:,tt)))
                            count_CA3_CA1_old(1,tt)=count_CA3_CA1_old(1,tt)+1;
                            phase_CA3_CA1_old_all{tt}(count_CA3_CA1_old(1,tt))=circ_mean(data{i}.coherence{CA1_ref(1)}.neuron(m).phi_old(:,tt));
                        end
                    end
                    if isempty(find(isnan(data{i}.coherence{CA1_ref(1)}.neuron(m).phi_old))) && isempty(find(isnan(data{i}.coherence{CA1_ref(1)}.neuron(m).phi_new)))
                        if circ_rtest(data{i}.coherence{CA1_ref(1)}.neuron(m).phi_old(:))<thr
                           if circ_rtest(data{i}.coherence{CA1_ref(1)}.neuron(m).phi_new(:))<thr
                            kk=kk+1;
                            phase_CA3_CA1(kk,:)=[circ_mean((data{i}.coherence{CA1_ref(1)}.neuron(m).phi_old)) circ_mean((data{i}.coherence{CA1_ref(1)}.neuron(m).phi_new))];
                        end
                    end
                    end
                end
    end
end

% run SVM
label(1:10)=0;
label(11:20)=1;
%%
% meanfrM=[phase_CA1_CA1'];
meanfrM=[phase_CA3_CA1'];
%meanfrM=[phase_CA3_CA3'];
trial=1:20
for i=1:length(trial)
    label_train=label';
    label_train(i)=[];
    label_train=squeeze(label_train);
    meanfr_train=meanfrM;
    meanfr_train(i,:)=[];
    meanfr_train=squeeze(meanfr_train);
     model = svmtrain(label_train, meanfr_train, params);
    [predicted_label , accuracy, decision_values] = svmpredict(label(i), squeeze(meanfrM(i,:)), model);
    test2(i)=predicted_label;
end
perf_all=length(find(test2==label))/length(trial)
%%
phase_CA3_CA1_old=phase_CA3_CA1(:,1:10);
phase_CA3_CA1_new=phase_CA3_CA1(:,11:20);
phase_CA1_CA1_new=phase_CA1_CA1(:,11:20);
phase_CA1_CA1_old=phase_CA1_CA1(:,1:10);
for i=1:size(phase_CA3_CA1,1)
    [h(i) p(i)]=ttest2(phase_CA3_CA1_old(i,:),phase_CA3_CA1_new(i,:));
end
phase_CA1_CA1_new=phase_CA1_CA1(:,11:20);
phase_CA1_CA1_old=phase_CA1_CA1(:,1:10);
for i=1:size(phase_CA1_CA1,1)
    [h_CA1(i) p_CA1(i)]=ttest2(phase_CA1_CA1_old(i,:),phase_CA1_CA1_new(i,:));
end



