function [result] = UW_spike_SWR_above_mean(filename,threshold,start_time,end_time)
    filter_swr = [150 250];

    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(filename,[1 1 1 1 1], 1, 1, [] );

    for i=1:length(Timestamps)
        for j=1:512
            TSArray(j,i) = Timestamps(i)+(1e6/SampleFrequencies(1,1)*(j-1));
        end
    end

    TSArray = TSArray(:)/1e6;
    csc_samples(:,1) = TSArray;
    csc_samples(:,2) = (Samples(:).*0.000000091552734375000002)*1e6;

    csc_samples = csc_samples(csc_samples(:,1)>=start_time,:);
    csc_samples = csc_samples(csc_samples(:,1)<=end_time,:);

    Fs = SampleFrequencies(1,1);
    n=2;
    Wn=[filter_swr(1) filter_swr(2)];
    Fn=Fs/2;
    ftype='bandpass';
    [b, a] = butter(n,Wn/Fn,ftype);

    csc_samples(:,3) = filtfilt(b,a,csc_samples(:,2));
    csc_samples(:,4) = abs(hilbert(csc_samples(:,3)));
    csc_samples(:,5) = smoothdata(csc_samples(:,4),'gaussian',40);
    
    SD3_gaussian = (std(csc_samples(:,5))*threshold) + mean(csc_samples(:,5));
    SD2_gaussian = (std(csc_samples(:,5))*2) + mean(csc_samples(:,5));

    flag=0;
    swr_num=1;

    for i=1:size(csc_samples,1)
        if(flag==0 & csc_samples(i,5)>=SD3_gaussian)
            swr_start = [i csc_samples(i,1)];
            flag = 1;
        elseif(flag==1 & csc_samples(i,5)<SD2_gaussian)
            swr_end = [i csc_samples(i,1)];
            swr(swr_num,:) = [swr_start swr_end];
            flag = 0;
            swr_num = swr_num+1;
        end
    end

    swr(:,5) = swr(:,4)-swr(:,2);
    swr = swr(swr(:,5)>=0.014,:);
    
    result.num_swr = size(swr,1);
    result.length_swr = mean(swr(:,5));
    result.swr = swr;
    result.csc_samples = csc_samples;
    result.sd_threshold = SD3_gaussian;
    result.mean = mean(csc_samples(:,5));
    result.filename = filename;
end