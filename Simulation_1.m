%% Exp0 Shell
    % Run Exp0 for several ISIs
    
    clear
% user-defined parameters ======> edit only this part
ISIs = .15:.05:.5;
for i = 1:length(ISIs)
    
    CSUS = ISIs(i);
    CS_duration = CSUS + .02;
    CS_rate = 100;
    USUS = 10;
    US_duration = .02;
    US_rate = 500;
    num_trials = 600;

% create parameters for run_experiment (everything is in seconds)
    trial_length                    = USUS; 
    CS_onset                        = 0.001;
    CS_offset                       = CS_duration;
    US_onset                        = CSUS;
    US_offset                       = CSUS+US_duration;
    archive0                        = zeros(1,800); 
    halt_at_pause                   = 0;
    
    if i == 1
        probe_start = 300;
        probe_freq = 20;
    else
        probe_start = Inf;
        probe_freq = 0;
        
    end
%% RUN : Acquisition

% run experiment
     [V{i},exc, inh, reserve, archive_tracker, record_activation_energy, readout_refractory_timer, readout_activation_energy,pause_trial] = ...
                            ...               
                            Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                                           US_onset, US_offset, US_rate,...
                                           trial_length, 0.5, num_trials,...
                                           archive0, halt_at_pause, probe_start, probe_freq);
                                       
    if i ==1
        V_spike = 10;
        record_activation_energy0  = 2.5;
        readout_activation_energy0 = 2;

         % electrophysiology on spikes
         figure() 
         spike_times = cell(num_trials,1);
         for n = 0:num_trials-1
           trial_voltage = V{i}(n+1,:);
           spike_indices = find(trial_voltage == V_spike);
           spike_times{n+1} = spike_indices/(1000);
           trial_height = (n+1)*ones(1,length(spike_indices));  
           if n+1 >= probe_start && mod(n+1,probe_freq) == 0
                   scatter(spike_times{n+1},trial_height,10,'ro','filled');  
           else   scatter(spike_times{n+1},trial_height,10,'ko','filled');  
           end
           hold on;  
         end

         for i = 1:length(CS_onset)
             line([CS_onset(i),CS_onset(i)],[0,num_trials],'Color','g','LineWidth',4)
             line([CS_offset(i),CS_offset(i)],[0,num_trials],'Color','g','LineWidth',4)
         end
         for j = 1:length(US_onset)
             line([US_onset(j),US_onset(j)],[0,num_trials],'Color','r','LineWidth',4)
         end
         hold off

         xlim([CS_onset(1),CS_offset(end)+0.5]);
         xlabel('Time (s)');
         ylabel('Trials');

%% RUN : Acquisition : Plot internal states
         if halt_at_pause
             whichtrial = pause_trial;
         else
             whichtrial = num_trials; 
         end

         %internal cell states
         figure() 
         subplot(2,1,1); plot(V{i}(whichtrial,:));
         subplot(2,1,2); plot(inh(whichtrial,:));
         

%%
        bin_width = 20;
        num_bins = floor(size(V{i},2)/bin_width);
        spikes = double(V{i}==10);
        collapse_spikes = sum(spikes,1);

        rate = zeros(1,num_bins);
        for j = 1:num_bins
            rate(j) = (1/num_trials)*(1/(bin_width*.001))*sum(collapse_spikes((j-1)*bin_width + 1: j*bin_width));
        end
        rate = smooth(rate); 

        figure()
        bar(bin_width*(1:num_bins),rate)
        xlabel('Time (ms)')
        ylabel('Rate (Hz)')
        xlim([0,1000])

    end
    pause
end 
%%        
for i = 1:length(V)
    
    experiment = V{i};
    bin_width = 20;
    num_bins = floor(size(experiment,2)/bin_width);
    spikes = double(experiment==10);
    collapse_spikes = sum(spikes,1);

    rate = zeros(1,num_bins);
    for j = 1:num_bins
        rate(j) = (1/num_trials)*(1/(bin_width*.001))*sum(collapse_spikes((j-1)*bin_width + 1: j*bin_width));
    end
    rate = smooth(rate); 

    
    [rate_max, max_loc] = max(rate);
    
    downward = rate(max_loc:end);
    
    subthreshold = find(downward<= 60) + max_loc - 1;
    
    pause_onset(i) = subthreshold(1); % pause onset
    
    [rate_min, min_loc(i)] = min(rate); % pause max
    
    pause_offset(i) = subthreshold(end); % pause_offset
    
    % Convert to correct timescale
    
    pause_onset(i) = bin_width*(pause_onset(i) - .5);
    min_loc(i) = bin_width*(min_loc(i) - .5);
    pause_offset(i) = bin_width*(pause_offset(i) - .5);
    
   figure;
   bar(bin_width*(1:num_bins), rate);
   xlim([0,1000])
   hold on;
   line([0 1000*USUS], [rate_max rate_max], 'Color', 'r', 'LineWidth',4)
   line([pause_onset(i) pause_onset(i)] + bin_width,[0,rate_max],'Color', 'g', 'LineWidth',4)
   line([pause_offset(i) pause_offset(i)] + bin_width,[0,rate_max],'Color', 'g', 'LineWidth',4)
   pause
   close;


end
  figure
  plot(ISIs*1000,min_loc,'o');
  hold on;
  plot(ISIs*1000, ISIs*1000, 'Color','r')
  xlabel('ISI (ms)')
  ylabel('Pause Maximum (ms)')
  
  figure
  plot(ISIs*1000,pause_onset,'o')
  xlabel('ISI (ms)')
  ylabel('Pause Onset (ms)')
  
  figure
  plot(ISIs*1000,pause_offset,'o')
  xlabel('ISI (ms)');
  ylabel('Pause Offset (ms)')
  
  figure
  plot(ISIs*1000,(pause_offset-pause_onset), 'o')
  xlabel('ISI (ms)');
  ylabel('Pause duration (ms)')

  
%% Test for CS duration effect:

%Experimental parameters
CSUS = .15;
CS_duration = CSUS + .1;
CS_rate = 100;
USUS = 10;
US_duration = .02;
US_rate = 500;
num_trials = 600;
trial_length                    = USUS; 
CS_onset                        = 0.001;
CS_offset                       = CS_duration;
US_onset                        = CSUS;
US_offset                       = CSUS+US_duration;
archive0                        = zeros(1,800); 
halt_at_pause                   = 0;

probe_start = Inf;
probe_freq = 0;

     [test_V,exc, inh, reserve, archive_tracker, record_activation_energy, readout_refractory_timer, readout_activation_energy,pause_trial] = ...
                            ...               
                            Purkinje_Cell(CS_onset, CS_offset, CS_rate,...
                                           US_onset, US_offset, US_rate,...
                                           trial_length, 0.5, num_trials,...
                                           archive0, halt_at_pause, probe_start, probe_freq);

% Plot

%% RUN : Acquisition : plot raster
V_spike = 10;

% electrophysiology on spikes
 figure() 
 spike_times = cell(num_trials,1);
 for n = 0:num_trials-1
   trial_voltage = V(n+1,:);
   spike_indices = find(trial_voltage == V_spike);
   spike_times{n+1} = spike_indices/(1000);
   trial_height = (n+1)*ones(1,length(spike_indices));  
   if n+1 >= probe_start && mod(n+1,probe_freq) == 0
           scatter(spike_times{n+1},trial_height,10,'ro','filled');  
   else   scatter(spike_times{n+1},trial_height,10,'ko','filled');  
   end
   hold on;  
 end
 
 for i = 1:length(CS_onset)
     line([CS_onset(i),CS_onset(i)],[0,num_trials],'Color','g','LineWidth',4)
     line([CS_offset(i),CS_offset(i)],[0,num_trials],'Color','g','LineWidth',4)
 end
 for j = 1:length(US_onset)
     line([US_onset(j),US_onset(j)],[0,num_trials],'Color','r','LineWidth',4)
 end
 hold off
 
 xlim([CS_onset(1),CS_offset(end)+0.5]);
 xlabel('Time (s)');
 ylabel('Trials');
