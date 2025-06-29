clear all; clc; close all;

for i=1:height(info)
    nlx(i).evt_filename = info.EVENT{i};
    nlx(i).hpc_filename = info.HIPP_NCS{i};
    nlx(i).pfc_filename = info.PFC_NTT{i};

    nlx(i).pre.time = [info.start_time1(i) info.end_time1(i)];
    nlx(i).pre.duration = info.end_time1(i) - info.start_time1(i);
    nlx(i).robot.time = [info.start_time2(i) info.end_time2(i)];
    nlx(i).robot.duration = info.end_time2(i) - info.start_time2(i);
    nlx(i).nest.time = [info.end_time2(i) info.end_time2(i)+600];
    nlx(i).nest.duration = nlx(i).nest.time(2) - nlx(i).nest.time(1);

    [~, ~, raw_data] = xlsread(nlx(i).evt_filename);
    nlx(i).evt = cell2mat(raw_data);
    nlx(i).evt(:,2) = nlx(i).evt(:,2)/1e6;

    pre_idx = nlx(i).evt(:,2) >= nlx(i).pre.time(1) & nlx(i).evt(:,2) <= nlx(i).pre.time(2);
    nlx(i).pre.evt = nlx(i).evt(pre_idx,:);

    robot_idx = nlx(i).evt(:,2) >= nlx(i).robot.time(1) & nlx(i).evt(:,2) <= nlx(i).robot.time(2);
    nlx(i).robot.evt = nlx(i).evt(robot_idx,:);

    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = ...
        Nlx2MatCSC(nlx(i).hpc_filename, [1 1 1 1 1], 1, 1, [] );

    for j = 1:length(Timestamps)
        for k = 1:512
            TSArray(k,j) = Timestamps(j) + (1e6/SampleFrequencies(1,1)*(k-1));
        end
    end
    TSArray = TSArray(:)/1e6;
    
    csc_samples(:,1) = TSArray;
    csc_samples(:,2) = Samples(:);
    
    [b, a] = butter(n, Wn/Fn, ftype);
    csc_samples(:,3) = filtfilt(b, a, csc_samples(:,2));
    hilbert_signal = hilbert(csc_samples(:,3));
    csc_samples(:,4) = angle(hilbert_signal);

    nlx(i).hpc.ncs_samples = csc_samples;
    clear Timestamps ChannelNumbers SampleFrequencies NumberOfValidSamples Samples Header TSArray csc_samples        

    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ... 
        Nlx2MatSpike(info.PFC_NTT{i}, [1 1 1 1 1], 1, 1, [] );

    temp = permute(Samples,[1 3 2]);
    nlx(i).pfc.spike_data = [Timestamps/1e6; CellNumbers; Features(1:4,:); temp(:,:,1); temp(:,:,2);temp(:,:,3);temp(:,:,4);]';
    nlx(i).pfc.spike_data = nlx(i).pfc.spike_data(nlx(i).pfc.spike_data(:,2)==info.Num_cell(i),:);

    clear Timestamps ScNumbers CellNumbers Features Samples Header temp

    spike_phases = zeros(size(nlx(i).pfc.spike_data, 1), 1);
    hpc_times = nlx(i).hpc.ncs_samples(:,1);
    hpc_phases = nlx(i).hpc.ncs_samples(:,4);
    
    for spike_idx = 1:size(nlx(i).pfc.spike_data, 1)
        spike_time = nlx(i).pfc.spike_data(spike_idx,1);
        [~, closest_idx] = min(abs(hpc_times - spike_time));
        
        if hpc_times(closest_idx) > spike_time && closest_idx > 1
            idx1 = closest_idx - 1;
            idx2 = closest_idx;
        elseif closest_idx < length(hpc_times)
            idx1 = closest_idx;
            idx2 = closest_idx + 1;
        else
            idx1 = closest_idx - 1;
            idx2 = closest_idx;
        end
        
        t1 = hpc_times(idx1);
        t2 = hpc_times(idx2);
        p1 = hpc_phases(idx1);
        p2 = hpc_phases(idx2);
        
        if abs(p2 - p1) > pi
            if p2 > p1
                p1 = p1 + 2*pi;
            else
                p2 = p2 + 2*pi;
            end
        end
        
        spike_phases(spike_idx) = p1 + (spike_time - t1) * (p2 - p1) / (t2 - t1);
        spike_phases(spike_idx) = mod(spike_phases(spike_idx) + pi, 2*pi) - pi;
    end
    
    nlx(i).pfc.spike_phases = [nlx(i).pfc.spike_data(:,1) spike_phases];
    clear spike_phases hpc_times hpc_phases

    pre_evt_phases = [];
    for evt_idx = 1:size(nlx(i).pre.evt, 1)
        evt_time = nlx(i).pre.evt(evt_idx,2);
        evt_start = evt_time - time_windows_before;
        evt_end = evt_time + time_windows_after;
        
        spike_times = nlx(i).pfc.spike_phases(:,1);
        spike_phase_values = nlx(i).pfc.spike_phases(:,2);
        window_spikes = spike_times >= evt_start & spike_times <= evt_end;
        pre_evt_phases{evt_idx} = spike_phase_values(window_spikes);
    end
    
    nlx(i).pre.spike_phases = pre_evt_phases;
    clear pre_evt_phases evt_time evt_start evt_end window_spikes

    robot_evt_phases = [];
    for evt_idx = 1:size(nlx(i).robot.evt, 1)
        evt_time = nlx(i).robot.evt(evt_idx,2);
        evt_start = evt_time - time_windows_before;
        evt_end = evt_time + time_windows_after;
        
        spike_times = nlx(i).pfc.spike_phases(:,1);
        spike_phase_values = nlx(i).pfc.spike_phases(:,2);
        window_spikes = spike_times >= evt_start & spike_times <= evt_end;
        robot_evt_phases{evt_idx} = spike_phase_values(window_spikes);
    end
    
    nlx(i).robot.spike_phases = robot_evt_phases;
    clear robot_evt_phases evt_time evt_start evt_end window_spikes

    if ~isempty(nlx(i).pre.spike_phases) && ~isempty(nlx(i).robot.spike_phases)
        combined_pre_phases{i} = cell2mat(nlx(i).pre.spike_phases');
        combined_robot_phases{i} = cell2mat(nlx(i).robot.spike_phases');
    end

    complex_phases_pre = exp(1j * combined_pre_phases{i});
    nlx(i).pre.plv = abs(mean(complex_phases_pre));
    nlx(i).pre.mean_phase = angle(mean(complex_phases_pre));

    complex_phases_robot = exp(1i * combined_robot_phases{i});
    nlx(i).robot.plv = abs(mean(complex_phases_robot));
    nlx(i).robot.mean_phase = angle(mean(complex_phases_robot));
    
    result.pfc_filename(i) = string(nlx(i).pfc_filename);
    result.hpc_filename(i) = string(nlx(i).hpc_filename);
    result.evt_filename(i) = string(nlx(i).evt_filename);
    result.num_cell(i) = info.Num_cell(i);
    result.pre_plv(i) = nlx(i).pre.plv;
    result.pre_mean_phase(i) = nlx(i).pre.mean_phase;
    result.robot_plv(i) = nlx(i).robot.plv;
    result.robot_mean_phase(i) = nlx(i).robot.mean_phase;

    if ~isempty(combined_pre_phases{i})
        [pre_counts,~] = histcounts(combined_pre_phases{i}, edges, 'Normalization', 'probability');
        pre_hists(i,:) = pre_counts;
        [p_rayleigh, z_rayleigh] = circ_rtest(combined_pre_phases{i});
        result.pre_rayleigh_p(i) = p_rayleigh;
        result.pre_rayleigh_z(i) = z_rayleigh;
        mrl = circ_r(combined_pre_phases{i});
        result.pre_mrl(i) = mrl;
        mean_phase = circ_mean(combined_pre_phases{i});
        result.pre_mean_phase(i) = mean_phase;
        result.pre_num_spikes(i) = length(combined_pre_phases{i});
    else
        pre_hists(i,:) = zeros(1,36);
        result.pre_rayleigh_p(i) = NaN;
        result.pre_rayleigh_z(i) = NaN;
        result.pre_mrl(i) = NaN;
        result.pre_mean_phase(i) = NaN;
        result.pre_num_spikes(i) = 0;
    end

    if ~isempty(combined_robot_phases{i})
        [robot_counts,~] = histcounts(combined_robot_phases{i}, edges, 'Normalization', 'probability');
        robot_hists(i,:) = robot_counts;
        [p_rayleigh, z_rayleigh] = circ_rtest(combined_robot_phases{i});
        result.robot_rayleigh_p(i) = p_rayleigh;
        result.robot_rayleigh_z(i) = z_rayleigh;
        mrl = circ_r(combined_robot_phases{i});
        result.robot_mrl(i) = mrl;
        mean_phase = circ_mean(combined_robot_phases{i});
        result.robot_mean_phase(i) = mean_phase;
        result.robot_num_spikes(i) = length(combined_robot_phases{i});
    else
        robot_hists(i,:) = zeros(1,36);
        result.robot_rayleigh_p(i) = NaN;
        result.robot_rayleigh_z(i) = NaN;
        result.robot_mrl(i) = NaN;
        result.robot_mean_phase(i) = NaN;
        result.robot_num_spikes(i) = 0;
    end

    result.pre_num_spikes(i) = length(combined_pre_phases{i});
    result.robot_num_spikes(i) = length(combined_robot_phases{i});
end

sig_robot_idx = find([result.robot_rayleigh_p] <= 0.05 & [result.robot_num_spikes] >= 40);

sig_pre_phases = [];
for idx = sig_robot_idx'
    sig_pre_phases = [sig_pre_phases; vertcat(nlx(idx).pre.spike_phases{:})];
end

sig_robot_phases = [];
for idx = sig_robot_idx'
    sig_robot_phases = [sig_robot_phases; vertcat(nlx(idx).robot.spike_phases{:})];
end

valid_robot_phases = [];
for i = 1:height(result)
    if result.robot_num_spikes(i) >= 40
        valid_robot_phases = [valid_robot_phases; result.robot_mean_phase(i)];
    end
end

[~, sorted_idx] = sort(result.robot_rayleigh_p);
min_p_idx = sorted_idx(5);

spike_phases = nlx(min_p_idx).robot.spike_phases;
if iscell(spike_phases)
    min_p_phases = vertcat(spike_phases{:});
    if size(min_p_phases,2) >= 2
        min_p_phases = min_p_phases(:,2);
    end
else
    min_p_phases = spike_phases(:,2);
end

[p, z] = circ_rtest(min_p_phases);
mean_angle = circ_mean(min_p_phases);
r = circ_r(min_p_phases);
