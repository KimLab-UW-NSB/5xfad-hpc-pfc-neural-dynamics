clear all; clc; close all;

params.Fs = 1280;
params.fpass = [2 12];
params.pad = 0;
params.tapers = [5 9];
params.err = 0;
params.trialave = 0;

info = readtable('mpfc_hpc_info_spike_w_unique.xlsx');

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
    Wn = [2 100];
    Fn = Fs/2;
    ftype = 'bandpass';
    [b, a] = butter(n, Wn/Fn, ftype);
    csc_samples(:,3) = filtfilt(b, a, csc_samples(:,2));

    nlx(i).hpc.ncs_samples = csc_samples;
    clear Timestamps ChannelNumbers SampleFrequencies NumberOfValidSamples Samples Header TSArray csc_samples

    nlx(i).pre.samples = nlx(i).hpc.ncs_samples(nlx(i).hpc.ncs_samples(:,1)>=info.start_time1(i)&nlx(i).hpc.ncs_samples(:,1)<=info.end_time1(i),3);
    nlx(i).robot.samples = nlx(i).hpc.ncs_samples(nlx(i).hpc.ncs_samples(:,1)>=info.start_time2(i)&nlx(i).hpc.ncs_samples(:,1)<=info.end_time2(i),3);

    [nlx(i).pre.S,nlx(i).pre.f] = mtspectrumc(nlx(i).pre.samples,params);
    [nlx(i).robot.S,nlx(i).robot.f] = mtspectrumc(nlx(i).robot.samples,params);
end

f_common = linspace(2, 12, 100);

pre_S_mean = zeros(size(f_common));
pre_S_table = array2table(f_common', 'VariableNames', {'Frequency'});
pre_S_percent_table = array2table(f_common', 'VariableNames', {'Frequency'});

for i = 1:length(nlx)
    S_interp = interp1(nlx(i).pre.f, nlx(i).pre.S, f_common, 'spline');
    window_size = 5;
    S_interp = movmean(S_interp, window_size);
    
    pre_S_mean = pre_S_mean + S_interp;
    pre_S_percent = (S_interp ./ sum(S_interp)) * 100;
    
    pre_S_table = addvars(pre_S_table, S_interp', 'NewVariableNames', {nlx(i).hpc_filename});
    pre_S_percent_table = addvars(pre_S_percent_table, pre_S_percent', ...
        'NewVariableNames', {[nlx(i).hpc_filename '_percent']});
        
    if i == 1
        pre_S_percent_mean = pre_S_percent;
    else
        pre_S_percent_mean = pre_S_percent_mean + pre_S_percent;
    end
end

pre_S_mean = pre_S_mean / length(nlx);
pre_S_percent_mean = pre_S_percent_mean / length(nlx);

robot_S_mean = zeros(size(f_common));
robot_S_table = array2table(f_common', 'VariableNames', {'Frequency'});
robot_S_percent_table = array2table(f_common', 'VariableNames', {'Frequency'});

for i = 1:length(nlx)
    S_interp = interp1(nlx(i).robot.f, nlx(i).robot.S, f_common, 'spline');
    window_size = 5;
    S_interp = movmean(S_interp, window_size);
    
    robot_S_mean = robot_S_mean + S_interp;
    robot_S_percent = (S_interp ./ sum(S_interp)) * 100;
    
    robot_S_table = addvars(robot_S_table, S_interp', 'NewVariableNames', {nlx(i).hpc_filename});
    robot_S_percent_table = addvars(robot_S_percent_table, robot_S_percent', ...
        'NewVariableNames', {[nlx(i).hpc_filename '_percent']});
        
    if i == 1
        robot_S_percent_mean = robot_S_percent;
    else
        robot_S_percent_mean = robot_S_percent_mean + robot_S_percent;
    end
end

robot_S_mean = robot_S_mean / length(nlx);
robot_S_percent_mean = robot_S_percent_mean / length(nlx);



filename = 'spectrum_analysis_results_hipp_w_session.xlsx';

writetable(pre_S_table, filename, 'Sheet', 'pre');
writetable(robot_S_table, filename, 'Sheet', 'robot');
writetable(pre_S_percent_table, filename, 'Sheet', 'pre_percent');
writetable(robot_S_percent_table, filename, 'Sheet', 'robot_percent');

disp(['Spectrum analysis results saved to ' filename]);