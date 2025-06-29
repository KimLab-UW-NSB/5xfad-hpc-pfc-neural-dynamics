function [ position ] = cal_vt_any( filename, t_start, t_end )
    opts = delimitedTextImportOptions("NumVariables", 4);
    
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    opts.VariableNames = ["TimeStamp", "ExtractedX", "ExtractedY", "ExtractedAngle"];
    opts.VariableTypes = ["double", "double", "double", "double"];
    
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    data = readtable(filename, opts);
    position = [data.TimeStamp data.ExtractedY data.ExtractedX data.ExtractedAngle];
    position = position(position(:,3)>0,:);
    
    if nargin >= 3
        t_start = t_start*1e6;
        t_end = t_end*1e6;
        position = position(position(:,1)>=t_start&position(:,1)<=t_end,:);
    end
end