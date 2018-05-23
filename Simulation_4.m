%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 4: Multiple Pauses %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

clear; 
close all;

% user-defined parameters ======> edit only this part
CSUS = [0.2;0.5];
CS_duration = CSUS + .02;
CS_rate = 100;
USUS = 30;
US_duration = .02;
US_rate = 500;
numtrials = 800;
halt_at_pause = 0;
ms_per_step = 1;
probe_start = 300; % since which trial
probe_freq = 20; % every # trials
 
% release_profile = .5 works. But it's internally 0 in RunExperiment

% create parameters for run_experiment (everything is in seconds)                              
    trial_length                    = USUS; 
    CS_onset                        = [0.001;0.001];
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
 