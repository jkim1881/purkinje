%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experiment 4: CS Invariance after training %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic


clear
% user-defined parameters ======> edit only this part
CSUS = 0.2;
CS_duration = 0.32;
CS_rate = 100;
USUS = 15;
US_duration = .02;
US_rate = 500;
numtrials = 600;
probe_start = Inf; % since which trial
probe_freq = 0; % every # trials
change_CS = 1;

% create parameters for run_experiment (everything is in seconds)                        
    trial_length                    = USUS; 
    CS_onset                        = 0.001;
    CS_offset                       = CS_duration;
    US_onset                        = CSUS;
    US_offset                       = CSUS + US_duration;
    archive0                        = zeros(1,800);    
    halt_at_pause                   = 0;
    probe_start                     = 300;
    probe_freq                      = 20;
%% RUN : Acquisition

% run experiment
     [V,exc, inh, reserve, archive_tracker, record_activation_energy, readout_refractory_timer, readout_activation_energy,pause_trial] = ...
                            ...               
                            Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                                           US_onset, US_offset, US_rate,...
                                           trial_length, 0.5,numtrials,archive0,halt_at_pause,...
                                           probe_start, probe_freq);


%% RUN : Acquisition : plot raster
 
V_spike = 10;
record_activation_energy0  = 2.5;
readout_activation_energy0 = 2;

 % electrophysiology on spikes
 figure() 
 spike_times = cell(numtrials,1);
 for n = 1:numtrials
   trial_voltage = V(n,:);
   spike_indices = find(trial_voltage == V_spike);
   spike_times{n} = spike_indices/(1000);
   trial_height = (n)*ones(1,length(spike_indices));  
   if n >= probe_start && mod(n,probe_freq) == 0
           scatter(spike_times{n},trial_height,10,'ro','filled');  
   else   scatter(spike_times{n},trial_height,10,'ko','filled');  
   end
   hold on;  
 end
 
 for i = 1:length(CS_onset)
     line([CS_onset(i),CS_onset(i)],[0,numtrials],'Color','g')
     line([CS_offset(i),CS_offset(i)],[0,numtrials],'Color','g')
 end
 for j = 1:length(US_onset)
     line([US_onset(j),US_onset(j)],[0,numtrials],'Color','r')
 end
 hold off
 
 xlim([CS_onset(1),CS_offset(end)+0.5]);
 xlabel('Time (s)');
 ylabel('Trials');
 title('Model Purkinje Cell Raster. Probe trials in red');
 
 %% PROBE 
whichtrial = 300; % at which 'stage' of learning do we perform probe trial?
 
CS_probe_duration = [0.1 : 0.3 : 1];
CS_probe_rate = [50 100 200];

numdurations = length(CS_probe_duration);
numrates = length(CS_probe_rate);

archive0 = archive_tracker(whichtrial,:);
 
V_probe = cell(numdurations,numrates);
for iduration = 1:numdurations
    for irate = 1:numrates
             [V_probe{iduration,irate},~, ~, ~, ~, ~, ~, ~, ~] = ...
                            ...               
                            Purkinje_Cell(0.001, CS_probe_duration(iduration), CS_probe_rate(irate),...
                                           [], [], 0,...
                                           trial_length, 1.5, 1 ,archive0,0,...
                                           [], []);
    end
end
%%
% generate figs
addpath('mollymolly-x-Axis-Brackets');
figure()
for iduration = 1:numdurations
    for irate = 1:numrates
        subplot(numrates+1,numdurations,numdurations*(irate-1) + iduration);
        plot(1:1:1200,V_probe{iduration,irate}(1:1200));
        axis([1 1200 -150 50]);
        set(gca,'fontsize',18)
        %draw lines
        for i = 1:length(CS_onset)
            p1 = line([1,1],[-150, 50],'Color','g','LineWidth',3);
            p2 = line([CS_probe_duration(iduration)*1000,CS_probe_duration(iduration)*1000],[-150, 50],'Color','g','LineWidth',3);
        end
        for j = 1:length(US_onset)
            p3 = line([US_onset(j)*1000,US_onset(j)*1000],[-150, 50],'Color','r');
        end
        if irate == numrates && iduration == 1
            xlabel('Time (s)', 'FontSize', 20,'FontWeight','bold');
            ylabel('Membrane Potential (mV)', 'FontSize', 20,'FontWeight','bold');
        end
    end
end
leg = legend([p1 p3],{'CS duration','ISI'}); 
set(leg,'FontSize',20)

%% another gif

    plot(1:1:1200,V_probe{3,2}(1:1200));
        axis([1 1200 -150 50]);
        set(gca,'fontsize',18)
                %draw lines
        for i = 1:length(CS_onset)
            p1 = line([1,1],[-150, 50],'Color','g','LineWidth',3);
            p2 = line([CS_probe_duration(3)*1000,CS_probe_duration(3)*1000],[-150, 50],'Color','g','LineWidth',3);
        end
        for j = 1:length(US_onset)
            p3 = line([US_onset(j)*1000,US_onset(j)*1000],[-150, 50],'Color','r');
        end
        title(strcat('Long probe (700ms)'));
        
 %%
  figure();
  plot(-6:-1:-20);
  locations = [6,8,1,2,8,14,8,15,4,5];
  labels = {'one','two','three','four','five'};
  [hLeft, hBottom, hRight] = xAxisBrackets(locations,'',[],'ordered');
  set(hBottom,'Color','m')
 