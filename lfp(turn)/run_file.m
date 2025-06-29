clear all; clc; close all;

params.Fs = 1280;
params.fpass = [2 12];
params.pad = 0;
params.tapers = [3 5];
params.err = 0;
params.trialave = 1;

time_windows_before = 3;
time_windows_after = 0;

info = readtable('mpfc_hpc_w.xlsx');

for i = 1:height(info)
    nlx(i).evt_filename = info.EVENT{i};
    nlx(i).hpc_filename = info.HIPP_NCS{i};
    nlx(i).pfc_filename = info.PFC_NCS{i};

    nlx(i).pre.time = [info.start_time1(i) info.end_time1(i)];
    nlx(i).pre.duration = info.end_time1(i) - info.start_time1(i);
    nlx(i).robot.time = [info.start_time2(i) info.end_time2(i)];
    nlx(i).robot.duration = info.end_time2(i) - info.start_time2(i);
    nlx(i).nest.time = [info.end_time2(i) info.end_time2(i)+600];
    nlx(i).nest.duration = nlx(i).nest.time(2) - nlx(i).nest.time(1);

    [~, ~, raw_data] = xlsread(nlx(i).evt_filename);
    
    nlx(i).evt = cell2mat(raw_data);
    nlx(i).evt(:,2) = nlx(i).evt(:,2)/1e6;

    pre_idx =  nlx(i).evt(:,2) >= nlx(i).pre.time(1) & nlx(i).evt(:,2) <= nlx(i).pre.time(2);
    nlx(i).pre.evt = nlx(i).evt(pre_idx,:);
    
    robot_idx = nlx(i).evt(:,2) >= nlx(i).robot.time(1) & nlx(i).evt(:,2) <= nlx(i).robot.time(2);
    nlx(i).robot.evt = nlx(i).evt(robot_idx,:);

    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = ...
        Nlx2MatCSC(nlx(i).hpc_filename, [1 1 1 1 1], 1, 1, []);

    for j = 1:length(Timestamps)
        for k = 1:512
            TSArray(k,j) = Timestamps(j) + (1e6/SampleFrequencies(1,1)*(k-1));
        end
    end
    TSArray = TSArray(:)/1e6;
    
    csc_samples(:,1) = TSArray;
    csc_samples(:,2) = (Samples(:) * 0.000000091552734375000002) * 1e6;
    
    Fs = SampleFrequencies(1,1);
    n = 3;
    Wn = [2 50];
    Fn = Fs/2;
    ftype = 'bandpass';
    [b, a] = butter(n, Wn/Fn, ftype);
    csc_samples(:,3) = filtfilt(b, a, csc_samples(:,2));

    nlx(i).hpc.ncs_samples = csc_samples;
    clear Timestamps ChannelNumbers SampleFrequencies NumberOfValidSamples Samples Header TSArray csc_samples

    pre_evt_samples = [];
    for evt_idx = 1:size(nlx(i).pre.evt, 1)
        evt_time = nlx(i).pre.evt(evt_idx,2);
        [~, closest_idx] = min(abs(nlx(i).hpc.ncs_samples(:,1) - evt_time));
        sample_idx = (closest_idx - Fs*time_windows_before):(closest_idx + Fs*time_windows_after);
        pre_evt_samples(:,evt_idx) = nlx(i).hpc.ncs_samples(sample_idx,3);
    end
    nlx(i).pre.hpc_evt_samples = pre_evt_samples;
    clear pre_evt_samples evt_time sample_idx

    robot_evt_samples = [];
    for evt_idx = 1:size(nlx(i).robot.evt, 1)
        evt_time = nlx(i).robot.evt(evt_idx,2);
        [~, closest_idx] = min(abs(nlx(i).hpc.ncs_samples(:,1) - evt_time));
        sample_idx = (closest_idx - Fs*time_windows_before):(closest_idx + Fs*time_windows_after);
        robot_evt_samples(:,evt_idx) = nlx(i).hpc.ncs_samples(sample_idx,3);
    end
    nlx(i).robot.hpc_evt_samples = robot_evt_samples;
    clear robot_evt_samples evt_time sample_idx
 
    [nlx(i).pre.S,nlx(i).pre.f] = mtspectrumc(nlx(i).pre.hpc_evt_samples,params);
    [nlx(i).robot.S,nlx(i).robot.f] = mtspectrumc(nlx(i).robot.hpc_evt_samples,params);
end

pre_S_mean = zeros(size(nlx(1).pre.S));
pre_S_table = array2table(nlx(1).pre.f', 'VariableNames', {'Frequency'});
pre_S_percent_table = array2table(nlx(1).pre.f', 'VariableNames', {'Frequency'});

for i = 1:length(nlx)
    pre_S_mean = pre_S_mean + nlx(i).pre.S;
    pre_S_percent = (nlx(i).pre.S ./ sum(nlx(i).pre.S)) * 100;
    
    pre_S_table = addvars(pre_S_table, nlx(i).pre.S, 'NewVariableNames', {nlx(i).hpc_filename});
    pre_S_percent_table = addvars(pre_S_percent_table, pre_S_percent, ...
        'NewVariableNames', {[nlx(i).hpc_filename '_percent']});
    if i == 1
        pre_S_percent_mean = pre_S_percent;
    else
        pre_S_percent_mean = pre_S_percent_mean + pre_S_percent;
    end
end

pre_S_mean = pre_S_mean / length(nlx);
pre_S_percent_mean = pre_S_percent_mean / length(nlx);

robot_S_mean = zeros(size(nlx(1).robot.S));
robot_S_table = array2table(nlx(1).robot.f', 'VariableNames', {'Frequency'});
robot_S_percent_table = array2table(nlx(1).robot.f', 'VariableNames', {'Frequency'});

for i = 1:length(nlx)
    robot_S_mean = robot_S_mean + nlx(i).robot.S;
    robot_S_percent = (nlx(i).robot.S ./ sum(nlx(i).robot.S)) * 100;
    robot_S_table = addvars(robot_S_table, nlx(i).robot.S, 'NewVariableNames', {nlx(i).hpc_filename});
    robot_S_percent_table = addvars(robot_S_percent_table, robot_S_percent, ...
        'NewVariableNames', {[nlx(i).hpc_filename '_percent']});
    if i == 1
        robot_S_percent_mean = robot_S_percent;
    else
        robot_S_percent_mean = robot_S_percent_mean + robot_S_percent;
    end
end

robot_S_mean = robot_S_mean / length(nlx);
robot_S_percent_mean = robot_S_percent_mean / length(nlx);

filename = 'spectrum_analysis_results_hipp_w_before_3.xlsx';

writetable(pre_S_table, filename, 'Sheet', 'pre');
writetable(robot_S_table, filename, 'Sheet', 'robot');
writetable(pre_S_percent_table, filename, 'Sheet', 'pre_percent');
writetable(robot_S_percent_table, filename, 'Sheet', 'robot_percent');

disp(['Spectrum analysis results saved to ' filename]);