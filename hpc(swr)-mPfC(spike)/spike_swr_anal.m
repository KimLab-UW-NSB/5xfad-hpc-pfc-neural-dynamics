clear all; clc; close all;

info = readtable('info_spike_swr_h_ej.xlsx');
SD = 4;
PSTH_window = .5;
bar_window = .01;

for i=1:size(info,1)
    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = ... 
        Nlx2MatSpike(info.PFC_NTT{i}, [1 1 1 1 1], 1, 1, [] );

    temp = permute(Samples,[1 3 2]);
    spike_data = [Timestamps; CellNumbers; Features(1:4,:); temp(:,:,1); temp(:,:,2);temp(:,:,3);temp(:,:,4);]';
    spike_data = spike_data(spike_data(:,2)==info.Num_cell(i),:);

    hpc_swr_pre = UW_spike_SWR_above_mean(info.HIPP_NCS{i},SD,info.start_time1(i),info.end_time1(i));
    hpc_swr_robot = UW_spike_SWR_above_mean(info.HIPP_NCS{i},SD,info.start_time2(i),info.end_time2(i));
    hpc_swr_nest = UW_spike_SWR_nest_above_mean(info.HIPP_NCS{i},SD);

    result{i}.FR.pre_spikes = spike_data(spike_data(:,1)/1e6>=info.start_time1(i) & spike_data(:,1)/1e6 <= info.end_time1(i),1);
    result{i}.FR.robot_spikes = spike_data(spike_data(:,1)/1e6>=info.start_time2(i) & spike_data(:,1)/1e6 <= info.end_time2(i),1);
    result{i}.FR.nest_spikes = spike_data(spike_data(:,1)/1e6>=hpc_swr_nest.csc_samples(end,1)-600,1);

    result{i}.FR.pre_dur = (info.end_time1(i)-info.start_time1(i));
    result{i}.FR.robot_dur = (info.end_time2(i)-info.start_time2(i)); 
    result{i}.FR.nest_dur = 600;

    result{i}.FR.pre_hz = size(result{i}.FR.pre_spikes,1)/result{i}.FR.pre_dur;
    result{i}.FR.robot_hz = size(result{i}.FR.robot_spikes,1)/result{i}.FR.robot_dur;
    result{i}.FR.nest_hz = size(result{i}.FR.nest_spikes,1)/result{i}.FR.nest_dur;

    for j=1:size(hpc_swr_pre.swr,1)
        st = hpc_swr_pre.swr(j,2) - PSTH_window;
        en = hpc_swr_pre.swr(j,2) + PSTH_window;
        PSTH_spikes_pre{j} = spike_data(spike_data(:,1)/1e6>=st & spike_data(:,1)/1e6 <= en,1)/1e6;
        PSTH_spikes_pre{j} = PSTH_spikes_pre{j} -hpc_swr_pre.swr(j,2);
        [PSTH_pre(j,:) edge] = histcounts(PSTH_spikes_pre{j},-PSTH_window:bar_window:PSTH_window);
    end

    for j=1:size(hpc_swr_robot.swr,1)
        st = hpc_swr_robot.swr(j,2) - PSTH_window;
        en = hpc_swr_robot.swr(j,2) + PSTH_window;
        PSTH_spikes_robot{j} = spike_data(spike_data(:,1)/1e6>=st & spike_data(:,1)/1e6 <= en,1)/1e6;
        PSTH_spikes_robot{j} = PSTH_spikes_robot{j} -hpc_swr_robot.swr(j,2);
        [PSTH_robot(j,:) edge] = histcounts(PSTH_spikes_robot{j},-PSTH_window:bar_window:PSTH_window);
    end

    for j=1:size(hpc_swr_nest.swr,1)
        st = hpc_swr_nest.swr(j,2) - PSTH_window;
        en = hpc_swr_nest.swr(j,2) + PSTH_window;
        PSTH_spikes_nest{j} = spike_data(spike_data(:,1)/1e6>=st & spike_data(:,1)/1e6 <= en,1)/1e6;
        PSTH_spikes_nest{j} = PSTH_spikes_nest{j} -hpc_swr_nest.swr(j,2);
        [PSTH_nest(j,:) edge] = histcounts(PSTH_spikes_nest{j},-PSTH_window:bar_window:PSTH_window);
    end

    result{i}.spike_data = spike_data;
    result{i}.PSTH_pre = PSTH_pre;
    result{i}.PSTH_robot = PSTH_robot;
    result{i}.PSTH_nest = PSTH_nest;
    result{i}.edge = edge;
    result{i}.hpc_swr_pre = hpc_swr_pre;
    result{i}.hpc_swr_robot = hpc_swr_robot;
    result{i}.hpc_swr_nest = hpc_swr_nest;
    result{i}.PSTH_spikes_pre = PSTH_spikes_pre;
    result{i}.PSTH_spikes_robot = PSTH_spikes_robot;
    result{i}.PSTH_spikes_nest = PSTH_spikes_nest;

    result{i}.pre.PSTH_sum = sum(PSTH_pre);
    result{i}.pre.PSTH_mean = mean(PSTH_pre);
    result{i}.pre.PSTH_sum_Hz = sum(PSTH_pre)/bar_window;
    result{i}.pre.PSTH_mean_Hz = mean(PSTH_pre)/bar_window;
    
    result{i}.pre.PSTH_mean_z = zscore(mean(PSTH_pre)/bar_window);
   
    result{i}.robot.PSTH_sum = sum(PSTH_robot);
    result{i}.robot.PSTH_mean = mean(PSTH_robot);
    result{i}.robot.PSTH_sum_Hz = sum(PSTH_robot)/bar_window;
    result{i}.robot.PSTH_mean_Hz = mean(PSTH_robot)/bar_window;
    result{i}.robot.PSTH_mean_z = zscore(mean(PSTH_robot)/bar_window);
    
    result{i}.nest.PSTH_sum = sum(PSTH_nest);
    result{i}.nest.PSTH_mean = mean(PSTH_nest);
    result{i}.nest.PSTH_sum_Hz = sum(PSTH_nest)/bar_window;
    result{i}.nest.PSTH_mean_Hz = mean(PSTH_nest)/bar_window;
    result{i}.nest.PSTH_mean_z = zscore(mean(PSTH_nest)/bar_window);
   
end

pre_hz_values = cellfun(@(x) x.FR.pre_hz, result);
robot_hz_values = cellfun(@(x) x.FR.robot_hz, result);
nest_hz_values = cellfun(@(x) x.FR.nest_hz, result);

pre_robot_corr = corr(pre_hz_values', robot_hz_values');
pre_nest_corr = corr(pre_hz_values', nest_hz_values');
nest_robot_corr = corr(nest_hz_values', robot_hz_values');