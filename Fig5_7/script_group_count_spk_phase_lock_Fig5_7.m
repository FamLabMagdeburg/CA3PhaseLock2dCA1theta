% Count how many neurons are phase-locked to local or CA1 theta and whose theta
% phase can discriminate repeated from non-repeated odors.
% Used for plotting proximodistal paper fig.5-7 
% Author: Shihpi
% Date: 23/Apr/2021
% Plotted the pie chart counting the proporiton of neurons which are phase
% locked, phase diff, phase locked but not diff, etc.
% Also plot the proportion of phase -locked neurons accross CA3
% proximodistal axis, although there is no clear preferential distribution.
% For fig. 7 use line 116 instead of line 115 for CA1 local phase locking counting and line 183 instead of line 182 for CA3. 


clear all
close all


cd /LFP_spectra
load region
count_CA1_spk_CA1_theta=zeros(2,5);
count_CA3_spk_CA3_theta=zeros(2,5);
CA1_count=0;
CA1_seg=[];
CA3_count=0;
CA3_seg=[];

thr=0.05;
aa=zeros(5,1);
bb=zeros(5,1);
res_CA1_spk_CA1_theta={};
for i=1:5
    for j=1:5
        res_CA1_spk_CA1_theta(i,j).mean_phase=0;
        res_CA1_spk_CA1_theta(i,j).p_rtest=NaN;
        aa(i,j)=0;
    end
end

res_CA3_spk_CA3_theta={};
for i=1:5
    for j=1:5
        
        res_CA3_spk_CA3_theta(i).mean_phase=0;
        res_CA3_spk_CA3_theta(i).p_rtest=NaN;
        bb(i,j)=0;
    end
end
loadfile='theta';


for animal_ID=[1:5]
    animal_ID
    cc=0;
    
    if animal_ID==1
        load(['LE82_coherence_' loadfile '_all.mat'],'data')
    elseif animal_ID==2
        load(['LE83_coherence_' loadfile '_all.mat'],'data')
    elseif animal_ID==3
        load(['LE84_coherence_' loadfile '_all.mat'],'data')
    elseif animal_ID==4
        load(['LE87_coherence_' loadfile '_all.mat'],'data')
    elseif animal_ID==5
        load(['LE46_coherence_' loadfile '_all.mat'],'data')
    end
    
    prop=[];
    CA1_tet=[];
    CA3_tet=[];
    CA1_tet=find(region(animal_ID).area(:,1)==1);
    CA3_tet=find(region(animal_ID).area(:,1)==3);
    CA1_tet_LFP=[];
    CA3_tet_LFP=[];
    % use dCA1 as ref
    if find(region(animal_ID).area(CA1_tet,2)==1) 
        CA1_tet_LFP=intersect(CA1_tet,find(region(animal_ID).area(:,2)==1));
        elseif find(region(animal_ID).area(CA1_tet,2)==2)
        CA1_tet_LFP=find(region(animal_ID).area(CA1_tet,2)==2);
        elseif ~isempty(find(region(animal_ID).area(CA1_tet,2)==3))
        CA1_tet_LFP=find(region(animal_ID).area(CA1_tet,2)==3);
    end
    
    % use pCA1 as ref
%     if find(region(animal_ID).area(CA1_tet,2)==5) 
%         CA1_tet_LFP=intersect(CA1_tet,find(region(animal_ID).area(:,2)==5));
%         elseif find(region(animal_ID).area(CA1_tet,2)==4)
%         CA1_tet_LFP=find(region(animal_ID).area(CA1_tet,2)==4);
%         elseif ~isempty(find(region(animal_ID).area(CA1_tet,2)==3))
%         CA1_tet_LFP=find(region(animal_ID).area(CA1_tet,2)==3);
%     end

    
    for i=1:length(data)
        tet_spk=data{i}.tet_spk;
        
        if find(CA1_tet==tet_spk)
            LFP_tet_CA1= CA1_tet_LFP(1);
            %LFP_tet_CA1= tet_spk;
            for m=1:length(data{i}.coherence{LFP_tet_CA1}.neuron)
                if ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).C_old) | ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).C_new)
                CA1_count=CA1_count+1;
                for n=1:5
                    %for p=1:5
                    if region(animal_ID).area(tet_spk,2)==n
                        count_CA1_spk_CA1_theta(1,n)=count_CA1_spk_CA1_theta(1,n)+1;
                        CA1_seg(CA1_count)=n;
                        aa(n)=aa(n)+1;
                        if ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_bsl)
                            phase_mean_CA1(n).p_rtest(aa(n),1)=circ_rtest(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_bsl(:));
                            phase_mean_CA1(n).mean_phase(aa(n),1)=circ_mean(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_bsl(:));
                        else
                            phase_mean_CA1(n).p_rtest(aa(n),1)=NaN;
                            phase_mean_CA1(n).mean_phase(aa(n),1)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_study)
                            phase_mean_CA1(n).p_rtest(aa(n),2)=circ_rtest(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_study(:));
                            phase_mean_CA1(n).mean_phase(aa(n),2)=circ_mean(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_study(:));
                        else
                            phase_mean_CA1(n).p_rtest(aa(n),2)=NaN;
                            phase_mean_CA1(n).mean_phase(aa(n),2)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_old)
                            
                            phase_mean_CA1(n).p_rtest(aa(n),3)=circ_rtest(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_old(:));
                            phase_mean_CA1(n).mean_phase(aa(n),3)=circ_mean(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_old(:));
                        else
                            phase_mean_CA1(n).p_rtest(aa(n),3)=NaN;
                            phase_mean_CA1(n).mean_phase(aa(n),3)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_new)
                            phase_mean_CA1(n).p_rtest(aa(n),4)=circ_rtest(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_new(:));
                            phase_mean_CA1(n).mean_phase(aa(n),4)=circ_mean(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_new(:));
                        else
                            phase_mean_CA1(n).p_rtest(aa(n),4)=NaN;
                            phase_mean_CA1(n).mean_phase(aa(n),4)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_old) & ~isempty(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_new)
                            if circ_cmtest(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_old(:),data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_new(:))<thr
                                
                                
                                phase_mean_CA1(n).mean_phase(aa(n),5)=circ_mean(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_old(:));
                                
                                phase_mean_CA1(n).mean_phase(aa(n),6)=circ_mean(data{i}.coherence{LFP_tet_CA1}.neuron(m).phi_new(:));
                            else
                                phase_mean_CA1(n).mean_phase(aa(n),5) = NaN;
                                phase_mean_CA1(n).mean_phase(aa(n),6) = NaN;
                            end
                        else
                            phase_mean_CA1(n).mean_phase(aa(n),5) = NaN;
                            phase_mean_CA1(n).mean_phase(aa(n),6) = NaN;
                            
                        end
                        
                    end
                end
               end
            end
        end
        
        if find(CA3_tet==tet_spk)
            LFP_tet_CA3=CA1_tet_LFP(1);
            %LFP_tet_CA3= tet_spk;

            for m=1:length(data{i}.coherence{LFP_tet_CA3}.neuron)
                if ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).C_old) | ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).C_new)

                CA3_count=CA3_count+1;
                for n=1:5
                    
                    if region(animal_ID).area(tet_spk,2)==n
                        count_CA3_spk_CA3_theta(1,n)=count_CA3_spk_CA3_theta(1,n)+1;
                        CA3_seg(CA3_count)=n;
                        bb(n)=bb(n)+1;
                        if ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_bsl)
                            phase_mean_CA3(n).p_rtest(bb(n),1)=circ_rtest(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_bsl(:));
                            phase_mean_CA3(n).mean_phase(bb(n),1)=circ_mean(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_bsl(:));
                        else
                            phase_mean_CA3(n).p_rtest(bb(n),1)=NaN;
                            phase_mean_CA3(n).mean_phase(bb(n),1)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_study)
                            phase_mean_CA3(n).p_rtest(bb(n),2)=circ_rtest(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_study(:));
                            phase_mean_CA3(n).mean_phase(bb(n),2)=circ_mean(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_study(:));
                        else
                            phase_mean_CA3(n).p_rtest(bb(n),2)=NaN;
                            phase_mean_CA3(n).mean_phase(bb(n),2)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_old)
                            phase_mean_CA3(n).p_rtest(bb(n),3)=circ_rtest(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_old(:));
                            phase_mean_CA3(n).mean_phase(bb(n),3)=circ_mean(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_old(:));
                        else
                            phase_mean_CA3(n).p_rtest(bb(n),3)=NaN;
                            phase_mean_CA3(n).mean_phase(bb(n),3)=NaN;
                        end
                        if ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_new)
                            phase_mean_CA3(n).p_rtest(bb(n),4)=circ_rtest(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_new(:));
                            phase_mean_CA3(n).mean_phase(bb(n),4)=circ_mean(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_new(:));
                        else
                            phase_mean_CA3(n).p_rtest(bb(n),4)=NaN;
                            phase_mean_CA3(n).mean_phase(bb(n),4)=NaN;
                        end

                        if ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_old) & ~isempty(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_new)
                            if circ_cmtest(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_old(:),data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_new(:))<thr
                                
                                
                                phase_mean_CA3(n).mean_phase(bb(n),5)=circ_mean(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_old(:));
                                
                                phase_mean_CA3(n).mean_phase(bb(n),6)=circ_mean(data{i}.coherence{LFP_tet_CA3}.neuron(m).phi_new(:));
                            else
                                phase_mean_CA3(n).mean_phase(bb(n),5) = NaN;
                                phase_mean_CA3(n).mean_phase(bb(n),6) = NaN;
                            end
                        else
                            phase_mean_CA3(n).mean_phase(bb(n),5) = NaN;
                            phase_mean_CA3(n).mean_phase(bb(n),6) = NaN;
                            
                        end
                    end
                end
            end
            
        end
        
        end
        
    end
    clear data
    
end
phase_mean_CA1_pvalue=phase_mean_CA1(1).p_rtest;
phase_mean_CA1_phase=phase_mean_CA1(1).mean_phase;
phase_mean_CA3_pvalue=phase_mean_CA3(1).p_rtest;
phase_mean_CA3_phase=phase_mean_CA3(1).mean_phase;
%%
for i=2:5
    phase_mean_CA3_pvalue=[phase_mean_CA3_pvalue;phase_mean_CA3(i).p_rtest];
    phase_mean_CA3_phase=[phase_mean_CA3_phase;phase_mean_CA3(i).mean_phase];
end

CA1_old=intersect(find(phase_mean_CA1_pvalue(:,3)<thr),find(~isnan(phase_mean_CA1_pvalue(:,3))));
CA1_new=intersect(find(phase_mean_CA1_pvalue(:,4)<thr),find(~isnan(phase_mean_CA1_pvalue(:,4))));
CA1_ret=unique([CA1_old; CA1_new]);
CA1_dis=sum(~isnan(phase_mean_CA1_phase(:,5)));
CA1_diff=find(~isnan(phase_mean_CA1_phase(:,5)));
CA1_diff_seg=CA1_seg(CA1_diff);
CA1_ret_seg=CA1_seg(CA1_ret);
% calculate how many neurons are phase-locked, phase-diff, only phase
% locked old, only phase locked new, phase locked both old and new, both
% and diff, both but same
CA1_res=zeros(226,1);
CA1_res(CA1_ret,1)=1;
CA1_res(CA1_diff,2)=1;
CA1_res(CA1_old,3)=1;
CA1_res(CA1_new,4)=1;
CA1_res((find(CA1_res(:,3)+CA1_res(:,4)==2)),5)=1;
CA1_res((find(CA1_res(:,2)+CA1_res(:,5)==2)),6)=1;
CA1_res(find(CA1_res(:,2)==0 & CA1_res(:,5)==1),7)=1;
phase_diff_old_only=length(find(CA1_res(:,2)==1 & CA1_res(:,3)==1& CA1_res(:,5)==0));
phase_diff_new_only=length(find(CA1_res(:,2)==1 & CA1_res(:,4)==1 & CA1_res(:,5)==0));
phase_same_old_new_lock=find(CA1_res(:,5)-CA1_res(:,2)==1);

%%
CA3_old=intersect(find(phase_mean_CA3_pvalue(:,3)<thr),find(~isnan(phase_mean_CA3_pvalue(:,3))));
CA3_new=intersect(find(phase_mean_CA3_pvalue(:,4)<thr),find(~isnan(phase_mean_CA3_pvalue(:,4))));
CA3_ret=unique([CA3_old; CA3_new]);
CA3_dis=sum(~isnan(phase_mean_CA3_phase(:,5)));
CA3_diff=find(~isnan(phase_mean_CA3_phase(:,5)));
CA3_diff_seg=CA3_seg(CA3_diff);
CA3_ret_seg=CA3_seg(CA3_ret);
% calculate how many neurons are phase-locked, phase-diff, only phase
% locked old, only phase locked new, phase locked both old and new, both
% and diff, both but same
CA3_res=zeros(226,1);
CA3_res(CA3_ret,1)=1;
CA3_res(CA3_diff,2)=1;
CA3_res(CA3_old,3)=1;
CA3_res(CA3_new,4)=1;
CA3_res((find(CA3_res(:,3)+CA3_res(:,4)==2)),5)=1;
CA3_res((find(CA3_res(:,2)+CA3_res(:,5)==2)),6)=1;
CA3_res(find(CA3_res(:,2)==0 & CA3_res(:,5)==1),7)=1;
phase_diff_old_only=length(find(CA3_res(:,2)==1 & CA3_res(:,3)==1& CA3_res(:,5)==0));
phase_diff_new_only=length(find(CA3_res(:,2)==1 & CA3_res(:,4)==1 & CA3_res(:,5)==0));
phase_same_old_new_lock=find(CA3_res(:,5)-CA3_res(:,2)==1);
% plot proportion pie chart
figure
ax = gca(); 
pieData = [24 42 69]./135; 
explode=[0 0 1];
h = pie(ax, pieData, explode); 

% ax is the handle to the pie chart axes
% to use a pre-defined colormap (in this example, 'Spring')
% h is the output to pie()
ax.Colormap = winter(numel(h)/2);
labels = {'phase same both locked','others','phase diff'};
legend(labels,'Location','southoutside','Orientation','horizontal')
title('proportion all phase locked CA3')
%
figure
ax = gca(); 
pieData = [7 13 49 ]./69; 
h = pie(ax, pieData); 

% ax is the handle to the pie chart axes
% to use a pre-defined colormap (in this example, 'Spring')
% h is the output to pie()
ax.Colormap = summer(numel(h)/2);
labels = {'old only','new only','both'};
legend(labels,'Location','southoutside','Orientation','horizontal')
title('proportion phase diff CA3')
