function [ position ] = cal_vt( filename, t_start, t_end )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    [VT_Timestamps, X, Y, Angles, Header] = Nlx2MatVT(filename, [1 1 1 1 0 0], 1, 1, [] );

    position = [VT_Timestamps; Y; X; Angles];
    position = position';
    position = position(position(:,3)>0,:);
    
    if nargin >= 3
        t_start = t_start*1e6;
        t_end = t_end*1e6;

        position = position(position(:,1)>=t_start&position(:,1)<=t_end,:);
    end

end

