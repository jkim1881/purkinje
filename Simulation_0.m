%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 0: Basic Learning and Extinction %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

clear; 
close all;

% user-defined parameters ======> edit only this part
CSUS = 0.2;
CS_duration = CSUS + .02;
CS_rate = 100;
USUS = 15;
US_duration = .02;
US_rate = 500;
numtrials = 400;
halt_at_pause = 0;
ms_per_step = 1;
probe_start = 300; % since which trial
probe_freq = 20; % every # trials
 
% release_profile = .5 works. But it's internally 0 in RunExperiment

% create parameters for run_experiment (everything is in seconds)                              
    trial_length                    = USUS; 
    CS_onset                        = 0.001;
    CS_offset                       = CS_duration;
    US_onset                        = CSUS;
    US_offset                       = CSUS+US_duration;
    archive0                        = zeros(1,800); 
    
%% RUN : Acquisition

% run experiment
     [V,exc, inh, reserve, archive_tracker, record_activation_energy, readout_refractory_timer, readout_activation_energy,pause_trial] = ...
                            ...               
                            Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                                           US_onset, US_offset, US_rate,...
                                           trial_length, 0.5, numtrials,archive0, halt_at_pause,probe_start, probe_freq);


%% RUN : Acquisition : plot raster
V_spike = 10;
record_activation_energy0  = 2.5;
readout_activation_energy0 = 2;

 % electrophysiology on spikes
 figure() 
 spike_times = cell(numtrials,1);
 for n = 0:numtrials-1
   trial_voltage = V(n+1,:);
   spike_indices = find(trial_voltage == V_spike);
   spike_times{n+1} = spike_indices/(1000);
   trial_height = (n+1)*ones(1,length(spike_indices));  
   if n+1 >= probe_start && mod(n+1,probe_freq) == 0
           scatter(spike_times{n+1},trial_height,30,'ro','filled');  
   else   scatter(spike_times{n+1},trial_height,10,'ko','filled');  
   end
   hold on;  
 end
 set(gca,'fontsize',18)
 
 for i = 1:length(CS_onset)
     line([CS_onset(i),CS_onset(i)],[0,numtrials],'Color','g','LineWidth',4)
     line([CS_offset(i),CS_offset(i)],[0,numtrials],'Color','g','LineWidth',4)
 end
 for j = 1:length(US_onset)
     line([US_onset(j),US_onset(j)],[0,numtrials],'Color','r','LineWidth',4)
 end
 hold off
 
 xlim([CS_onset(1),CS_offset(end)+0.5]);
 xlabel('Time (s)', 'FontSize', 20,'FontWeight','bold');
 ylabel('Trials', 'FontSize', 20,'FontWeight','bold');
 title('Model Purkinje Cell Raster. Probe trials in red');
 
  %% RUN : Acquisition : Plot internal states
 if halt_at_pause
     whichtrial = pause_trial;
 else
     whichtrial = numtrials; 
 end
 
 %internal cell states
 figure() 
subplot(3,1,1); plot(0.001:0.001:size(V,2)/1000,archive_tracker(whichtrial,1:size(V,2)),'LineWidth',2);hold on
    line([CSUS,CSUS],[0,0.1],'LineStyle','--','Color','r','LineWidth',1);hold off
     set(gca,'fontsize',18)
     xlim([0 size(V,2)/1000])
    title('Archive', 'FontSize', 20,'FontWeight','bold');
     ylabel('Recorder Units', 'FontSize', 20,'FontWeight','bold');
subplot(3,1,2); plot(0.001:0.001:size(V,2)/1000,inh(whichtrial,1:size(V,2)),'LineWidth',2);hold on
    line([CSUS,CSUS],[0,15],'LineStyle','--','Color','r','LineWidth',1);hold off
     set(gca,'fontsize',18)
      xlim([0 size(V,2)/1000])
    title('Inhibition', 'FontSize', 20,'FontWeight','bold')
     ylabel('Inhibitory Current (mV/ms)', 'FontSize', 20,'FontWeight','bold');
subplot(3,1,3); plot(0.001:0.001:size(V,2)/1000,V(whichtrial,1:size(V,2)),'LineWidth',2);hold on
    line([CSUS,CSUS],[-150,50],'LineStyle','--','Color','r','LineWidth',1);hold off
     set(gca,'fontsize',18)
      xlim([0 size(V,2)/1000])
    title('Membrane Potential', 'FontSize', 20,'FontWeight','bold')
     ylabel('Trials', 'FontSize', 20,'FontWeight','bold');
     xlabel('Time (s)', 'FontSize', 20,'FontWeight','bold');

 
 
  %% RUN : Acquisition : Plot development of archive for different ISIs

whichtrial = [20 40 80 160 320 640]; 
isiss = [0.2:0.2:0.8];
archive0 = zeros(1,1500); 
 %internal cell states
 k= 0;
 inh_temp = cell(4,1);
 archive_tracker_temp = cell(4,1);
for isis = isiss
    k = k+1;
    % run for each isi
    [~,~, inh_temp{k}, ~, archive_tracker_temp{k}, ~, ~, ~,~] = ...
                            ...               
                            Purkinje_Cell(CS_onset, 0.2, 100,...
                                           isis, isis+US_duration, US_rate,...
                                           trial_length, 0.5, 800,archive0, 0,probe_start, probe_freq);
end


% plot figure
figure()
ylim = 0.14;
for k = 1:length(isiss)
    % plot
    subplot(length(isiss)+1,1,k);
    for itrial = 1:length(whichtrial)
        plot((1:1500)./1000, [archive_tracker_temp{k}(whichtrial(itrial),:) zeros(1,1500-size(archive_tracker_temp{k},2))],'LineWidth',2);hold on
        axis([0 1.5 0 ylim]);
         set(gca,'fontsize',18)
    end
    line([isiss(k),isiss(k)],[0,ylim],'LineStyle','--','Color','r','LineWidth',1);hold off
    if k == length(isiss)
        xlabel('Time (s)','FontSize', 20,'FontWeight','bold');
        ylabel('Recorder Units','FontSize', 20,'FontWeight','bold')
    end
end

%legend
labels = cell(1,length(whichtrial));
for itrial = 1:length(whichtrial)
    labels{itrial} = strcat('# trials = ',num2str(whichtrial(itrial)));
end
leg = legend(labels); 
set(leg,'FontSize',20)

 
%%
bin_width = 20;
num_bins = floor(size(V,2)/bin_width);
spikes = double(V==10);
collapse_spikes = sum(spikes,1);

rate = zeros(1,num_bins);
for j = 1:num_bins
    rate(j) = (1/numtrials)*(1/(bin_width*.001))*sum(collapse_spikes((j-1)*bin_width + 1: j*bin_width));
end
rate = smooth(rate); 

figure()
bar(bin_width*(1:num_bins),rate)
xlim([0,1000])
 
 
%  %% RUN : Extinction
% 
%   CSUS = 0.5;
%   CS_duration = 0.12;
%   CS_rate = 100;
%   USUS = 3;
%   US_duration = .02;
%   US_rate = 500;
%   numtrials = 90;
%  
%   % create parameters for run_experiment (everything is in seconds)                            
%       trial_length                    = USUS; 
%       CS_onset                        = 0.001;
%       CS_offset                       = CS_duration;
%       US_onset                        = 0;
%       US_offset                       = 0;
%       archive0                        = archive_tracker(end,:);  
%       
%    
%        [V,exc, inh, reserve, archive_tracker, record_activation_energy, readout_refractory_timer, readout_activation_energy,pause_trial] = ...
%                               ...               
%                              run_experiment(CS_onset, CS_offset, CS_rate,...
%                                            US_onset, US_offset, US_rate,...
%                                            trial_length, numtrials,archive0, halt_at_pause);
% 
% 
%   %% RUN : Extinction : plot raster
%   probe_start = Inf; % since which trial
%   probe_freq = Inf; % every # trials
%    
%   V_spike = 10;
%   record_activation_energy0  = 2.5;
%   readout_activation_energy0 = 2;
%   
%    % electrophysiology on spikes
%    figure() 
%    spike_times = cell(numtrials,1);
%    for n = 0:numtrials-1
%     trial_voltage = V(n*(trial_length*1000)+1 : (n+1)*(trial_length*1000));
%      spike_indices = find(trial_voltage == V_spike);
%      spike_times{n+1} = spike_indices/(1000);
%      trial_height = (n+1)*ones(1,length(spike_indices));  
%      if n+1 >= probe_start && mod(n+1,probe_freq) == 0
%              scatter(spike_times{n+1},trial_height,10,'ro','filled');  
%      else   scatter(spike_times{n+1},trial_height,10,'ko','filled');  
%      end
%      hold on;  
%    end
%    
%    for i = 1:length(CS_onset)
%        line([CS_onset(i),CS_onset(i)],[0,numtrials],'Color','g')
%        line([CS_offset(i),CS_offset(i)],[0,numtrials],'Color','g')
%    end
% 
%    line([.500,.500],[0,numtrials],'Color','r')
%    hold off
%    
%    xlim([CS_onset(1),CS_offset(end)+1]);
%    xlabel('Time (s)');
%    ylabel('Trials');
%    title('Model Purkinje Cell Raster. Probe trials in red');
%    
%    %% RUN : Extinction : Plot internal states
%    whichtrial = 13; 
%    
%    plotinterval = trial_length*(1000/ms_per_step)*(whichtrial-1)+...
%                  (CS_onset(1)*(1000/ms_per_step):1:(CS_offset(end)+1)*(1000/ms_per_step));
%    %internal cell states
%    figure() 
%    subplot(2,1,1); plot(plotinterval,V(plotinterval));title('pkj cell potential')
%    subplot(2,1,2); plot(plotinterval,inh(plotinterval));title('inhibitory current')
%    
%    %hypothetical states
%    figure() 
%    subplot(4,1,1); plot(plotinterval,reserve(plotinterval));title('recordable reserve')
%    subplot(4,1,2); plot(plotinterval,record_activation_energy(plotinterval));title('record E(activation) - threshold toggle');ylim([0 record_activation_energy0]);
%    subplot(4,1,3); plot(plotinterval,readout_refractory_timer(plotinterval));title('readout refractory')
%    subplot(4,1,4); plot(plotinterval,readout_activation_energy(plotinterval));title('readout E(activation) - threshold refractory');ylim([0 readout_activation_energy0]);
%    
%    
%    %% RUN : Extinction : PLOT archive tracker
%    figure()
%    for itrial = 1:numtrials
%        plot(1:size(archive_tracker,2),archive_tracker(itrial,:));ylim([0 0.1]);
%        title(['Trial# ',num2str(itrial)]);
%        pause(0.1)
%    end
%    