function [result] = UW_SWR_locomotion(filename_ncs, filename_nvt, threshold, start_time, end_time, session, risky_side, smooth_time)
    filter_low = [2 50];
    filter_swr = [150 250];

    pos = cal_vt(filename_nvt);
    pos(:,1) = pos(:,1)/1e6;
    
    if(strcmp(session,'pre'))
        start_time_pos = start_time;
        end_time_pos = end_time;
    elseif(strcmp(session,'robot'))
        start_time_pos = start_time;
        end_time_pos = end_time;
    else
        end_time_nest = pos(end,1);
        start_time_pos = end_time_nest - 600;
        if(start_time_pos <= end_time)
            start_time_pos = end_time;
            end_time_pos = end_time_nest;
        else
            end_time_pos = end_time_nest;
        end
    end

    pos = pos(pos(:,1)>=start_time_pos & pos(:,1)<=end_time_pos,:);

    x_pixel_size = 0.1508;
    y_pixel_size = 0.1371;
    
    pos(:,5) = 0;
    pos(:,6) = 0;
    pos(:,7) = 0;
    
    dt = diff(pos(:,1));
    dx = diff(pos(:,3));
    dy = diff(pos(:,2));
    
    distance_cm = sqrt(dx.^2 + dy.^2) * sqrt(x_pixel_size^2 + y_pixel_size^2);
    
    valid_dt = dt > 0;
    velocity = zeros(size(dt));
    velocity(valid_dt) = distance_cm(valid_dt) ./ dt(valid_dt);
    
    pos(2:end,5) = velocity;
    
    window_size = smooth_time * 30;
    smoothing_method = 'gaussian';
    pos(:,6) = smoothdata(pos(:,5), smoothing_method, window_size);
    
    velocity_threshold = 5;
    pos(:,7) = pos(:,6) > velocity_threshold;
    
    not_moving_periods = [];
    is_not_moving = false;
    start_idx = 0;
        
    min_not_moving_duration = 0.1;

    for i = 1:size(pos, 1)
        if pos(i, 6) <= velocity_threshold
            if ~is_not_moving
                start_idx = i;
                is_not_moving = true;
            end
        else
            if is_not_moving
                period_start_time = pos(start_idx, 1);
                period_end_time = pos(i-1, 1);
                
                if (period_end_time - period_start_time) >= min_not_moving_duration
                    not_moving_periods = [not_moving_periods; period_start_time, period_end_time];
                end

                is_not_moving = false;
            end
        end
    end
    
    if is_not_moving
        period_start_time = pos(start_idx, 1);
        period_end_time = pos(end, 1);
    
        if (period_end_time - period_start_time) >= min_not_moving_duration
            not_moving_periods = [not_moving_periods; period_start_time, period_end_time];
        end
    end
    
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = ...
        Nlx2MatCSC(filename_ncs, [1 1 1 1 1], 1, 1, [] );

    for i = 1:length(Timestamps)
        for j = 1:512
            TSArray(j,i) = Timestamps(i) + (1e6/SampleFrequencies(1,1)*(j-1));
        end
    end
    TSArray = TSArray(:)/1e6;
    
    csc_samples(:,1) = TSArray;
    csc_samples(:,2) = (Samples(:) * 0.000000091552734375000002) * 1e6;

    if(strcmp(session,'pre'))
    elseif(strcmp(session,'robot'))
    else
        end_time_nest = csc_samples(end,1);
        start_time = end_time_nest - 600;
        if(start_time <= end_time)
            start_time = end_time;
            end_time = end_time_nest;
        else
            end_time = end_time_nest;
        end
    end
    csc_samples = csc_samples(csc_samples(:,1) >= start_time & csc_samples(:,1) <= end_time, :);

    Fs = SampleFrequencies(1,1);
    n = 2;
    Wn = filter_swr;
    Fn = Fs/2;
    ftype = 'bandpass';
    [b, a] = butter(n, Wn/Fn, ftype);
    csc_samples(:,3) = filtfilt(b, a, csc_samples(:,2));

    csc_samples(:,4) = abs(hilbert(csc_samples(:,3)));
    csc_samples(:,5) = smoothdata(csc_samples(:,4), 'gaussian', 40);
    
    SD_threshold_gaussian = (std(csc_samples(:,5)) * threshold) + mean(csc_samples(:,5));
    SD2_gaussian = (std(csc_samples(:,5)) * 2) + mean(csc_samples(:,5));
    
    flag = 0;
    swr_num = 1;
    for i = 1:size(csc_samples,1)
        if flag == 0 && csc_samples(i,5) >= SD_threshold_gaussian
            swr_start = [i, csc_samples(i,1)];
            flag = 1;
        elseif flag == 1 && csc_samples(i,5) < SD2_gaussian
            swr_end = [i, csc_samples(i,1)];
            swr(swr_num,:) = [swr_start, swr_end];
            flag = 0;
            swr_num = swr_num + 1;
        end
    end
    swr(:,5) = swr(:,4) - swr(:,2);
    swr = swr(swr(:,5) >= 0.014, :);
    
    swr_pos_x_start = zeros(size(swr,1), 1);
    swr_pos_y_start = zeros(size(swr,1), 1);
    swr_pos_x_end = zeros(size(swr,1), 1);
    swr_pos_y_end = zeros(size(swr,1), 1);
    
    time_differences_start = abs(pos(:,1) - swr(:,2)');
    [~, closest_indices_start] = min(time_differences_start, [], 1);
    swr_pos_x_start = pos(closest_indices_start, 2);
    swr_pos_y_start = pos(closest_indices_start, 3);
    
    time_differences_end = abs(pos(:,1) - swr(:,4)');
    [~, closest_indices_end] = min(time_differences_end, [], 1);
    swr_pos_x_end = pos(closest_indices_end, 2);
    swr_pos_y_end = pos(closest_indices_end, 3);

    swr_with_pos = [swr, swr_pos_x_start, swr_pos_y_start, swr_pos_x_end, swr_pos_y_end];
    
    csc_samples_table = array2table(csc_samples, ...
        'VariableNames', {'Timestamp', 'RawSignal', 'FilteredSignal', 'HilbertAmplitude', 'SmoothSignal'});
    
    swr_table = array2table(swr_with_pos, ...
        'VariableNames', {'StartIndex', 'StartTime', 'EndIndex', 'EndTime', 'Duration', ...
        'StartPosX', 'StartPosY', 'EndPosX', 'EndPosY'});
    
    pos_table = array2table(pos, ...
        'VariableNames', {'Timestamp', 'X', 'Y', 'Angle', 'InstVelocity(cm/sec)', 'SmoothedVelocity(cm/sec)', 'MovementState'});
    
    if endsWith(filename_nvt, ".nvt", 'IgnoreCase', true)
        location = repmat({'undefined'}, height(swr_table), 1);
        
        idx = swr_table.StartPosY >= 560;
        location(idx) = {'nest'};
        
        idx = swr_table.StartPosY < 560 & swr_table.StartPosX >= 570 & swr_table.StartPosX <= 730;
        location(idx) = {'center'};
        
        idx = swr_table.StartPosY < 560 & swr_table.StartPosY <= 275 & swr_table.StartPosX < 570;
        location(idx) = {'left arm'};
        
        idx = swr_table.StartPosY < 560 & swr_table.StartPosY > 275 & swr_table.StartPosX < 570;
        location(idx) = {'left distal'};
        
        idx = swr_table.StartPosY < 560 & swr_table.StartPosY <= 275 & swr_table.StartPosX > 730;
        location(idx) = {'right arm'};
        
        idx = swr_table.StartPosY < 560 & swr_table.StartPosY > 275 & swr_table.StartPosX > 730;
        location(idx) = {'right distal'};
    elseif endsWith(filename_nvt, ".xlsx", 'IgnoreCase', true)
        location = repmat({'undefined'}, height(swr_table), 1);
        
        idx = swr_table.StartPosY >= 772;
        location(idx) = {'nest'};
        
        idx = swr_table.StartPosY < 772 & swr_table.StartPosX >= 847 & swr_table.StartPosX <= 1095;
        location(idx) = {'center'};
        
        idx = swr_table.StartPosY < 772 & swr_table.StartPosY <= 375 & swr_table.StartPosX < 847;
        location(idx) = {'left arm'};
        
        idx = swr_table.StartPosY < 772 & swr_table.StartPosY > 375 & swr_table.StartPosX < 847;
        location(idx) = {'left distal'};
        
        idx = swr_table.StartPosY < 772 & swr_table.StartPosY <= 375 & swr_table.StartPosX > 1095;
        location(idx) = {'right arm'};
        
        idx = swr_table.StartPosY < 772 & swr_table.StartPosY > 375 & swr_table.StartPosX > 1095;
        location(idx) = {'right distal'};
    end
    
    swr_table = addvars(swr_table, location, 'After', 'EndPosY', 'NewVariableNames', 'swr_location');
    
    movement_state = repmat({'moving'}, height(swr_table), 1);

    if exist('not_moving_periods', 'var') && ~isempty(not_moving_periods)
        for i = 1:height(swr_table)
            swr_start_time = double(swr_table.StartTime(i));
            
            for j = 1:size(not_moving_periods, 1)
                start_period = double(not_moving_periods(j, 1));
                end_period = double(not_moving_periods(j, 2));
                
                if swr_start_time >= start_period && swr_start_time <= end_period
                    movement_state{i} = 'not moving';
                    break;
                end
            end
        end
    end
    
    swr_table = addvars(swr_table, movement_state, 'After', 'swr_location', 'NewVariableNames', 'movement_state');
    
    if exist('risky_side', 'var') && ~isempty(risky_side)
        risky_side_column = repmat({risky_side}, height(swr_table), 1);
        swr_table = addvars(swr_table, risky_side_column, 'After', width(swr_table), 'NewVariableNames', 'risky_side');

        risky_location = repmat("undefined", height(swr_table), 1);

        if strcmp(risky_side, 'L')
            idx = ismember(location, {'left arm'}); risky_location(idx) = "risky arm";
            idx = ismember(location, {'left distal'}); risky_location(idx) = "risky distal";
            idx = ismember(location, {'right arm'}); risky_location(idx) = "safe arm";
            idx = ismember(location, {'right distal'}); risky_location(idx) = "safe distal";
            idx = ismember(location, {'nest'}); risky_location(idx) = "nest";
            idx = ismember(location, {'center'}); risky_location(idx) = "center";
        elseif strcmp(risky_side, 'R')
            idx = ismember(location, {'right arm'}); risky_location(idx) = "risky arm";
            idx = ismember(location, {'right distal'}); risky_location(idx) = "risky distal";
            idx = ismember(location, {'left arm'}); risky_location(idx) = "safe arm";
            idx = ismember(location, {'left distal'}); risky_location(idx) = "safe distal";
            idx = ismember(location, {'nest'}); risky_location(idx) = "nest";
            idx = ismember(location, {'center'}); risky_location(idx) = "center";
        end

        swr_table = addvars(swr_table, risky_location, 'After', width(swr_table), 'NewVariableNames', 'risky_location');
    end
    
    result.num_swr = size(swr,1);
    result.length_swr = mean(swr(:,5));
    result.swr = swr_table;
    result.csc_samples = csc_samples_table;
    result.sd_threshold = SD_threshold_gaussian;
    result.mean = mean(csc_samples(:,5));
    result.filename_ncs = filename_ncs;
    result.filename_nvt = filename_nvt;
    result.pos = pos_table;
    result.cm = array2table([x_pixel_size y_pixel_size], 'VariableNames', {'pixel_to_cm_X','pixel_to_cm_Y'});
    result.session = session;
end