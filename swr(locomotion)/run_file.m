clear all; clc; close all;

info = readtable('info_spike_swr_h3.xlsx');
threshold = 4;
smooth_time = 1;

for i=1:size(info,1)
    result{i}.pre = UW_SWR_locomotion(info.HIPP_NCS{i},info.NVT{i},threshold,info.start_time1(i),info.end_time1(i),'pre',info.risky_side{i},smooth_time);
    result{i}.robot = UW_SWR_locomotion(info.HIPP_NCS{i},info.NVT{i},threshold,info.start_time2(i),info.end_time2(i),'robot',info.risky_side{i},smooth_time);
    result{i}.nest = UW_SWR_locomotion(info.HIPP_NCS{i},info.NVT{i},threshold,0,info.end_time2(i),'nest',info.risky_side{i},smooth_time);

    writetable(result{i}.pre.swr,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','pre');
    writetable(result{i}.robot.swr,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','robot');
    writetable(result{i}.nest.swr,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','nest');
    writetable(result{i}.pre.summary_table,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','pre_summary');
    writetable(result{i}.robot.summary_table,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','robot_summary');
    writetable(result{i}.nest.summary_table,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','nest_summary');
end

clear all; clc; close all;

info = readtable('info_spike_swr_w3.xlsx');
threshold = 4;
smooth_time = 1;

for i=1:size(info,1)
    result{i}.pre = UW_SWR_locomotion(info.HIPP_NCS{i},info.NVT{i},threshold,info.start_time1(i),info.end_time1(i),'pre',info.risky_side{i},smooth_time);
    result{i}.robot = UW_SWR_locomotion(info.HIPP_NCS{i},info.NVT{i},threshold,info.start_time2(i),info.end_time2(i),'robot',info.risky_side{i},smooth_time);
    result{i}.nest = UW_SWR_locomotion(info.HIPP_NCS{i},info.NVT{i},threshold,0,info.end_time2(i),'nest',info.risky_side{i},smooth_time);

    writetable(result{i}.pre.swr,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','pre');
    writetable(result{i}.robot.swr,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','robot');
    writetable(result{i}.nest.swr,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','nest');
    writetable(result{i}.pre.summary_table,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','pre_summary');
    writetable(result{i}.robot.summary_table,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','robot_summary');
    writetable(result{i}.nest.summary_table,strcat(info.HIPP_NCS{i}(1:end-4),'.xlsx'),'Sheet','nest_summary');
end