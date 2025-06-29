clear all;clc; close all;

vt_name = 'H18_010422_VT1.nvt';

pos = cal_vt(vt_name);

%% Y-axis : Nest out_timing_return
fig = figure; %figure('menubar','none');

% % Add menus with Accelerators
% mymenu = uimenu('Parent',fig,'Label','Hot Keys');
% uimenu('Parent',mymenu,'Label','Zoom','Accelerator','z','Callback',@(src,evt)zoom(fig,'on'));
% uimenu('Parent',mymenu,'Label','Rotate','Accelerator','r','Callback',@(src,evt)rotate3d(fig,'on'));
% uimenu('Parent',mymenu,'Label','Pan','Accelerator','p','Callback',@(src,evt)pan(fig,'on'));

set(gcf,'units','normalized','outerposition',[0 0 1 1]);
hold on;
plot([700 700],[pos(1,1) pos(end,1)],'r--');hold on;
plot([600 600],[pos(1,1) pos(end,1)],'r--');
plot(pos(:,3),pos(:,1),'k.');

dcm = datacursormode(fig);
dcm.Enable = 'on';

disp('Click line to display a data tip, then press "Return"')
pause

vals = getCursorInfo(dcm);

for i=1:size(vals,2)
    temp(i,:) = [vals(i).Position vals(i).DataIndex];
end

nest_out = sortrows(temp,3,'ascend');


ref_time = 15 * 2;

figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);hold on;
plot(pos(:,3),pos(:,2),'.','color',([0.5 0.5 0.5]));
for i=1:size(nest_out,1)
    plot(pos(nest_out(i,3)-ref_time : nest_out(i,3)+ref_time,3),pos(nest_out(i,3)-ref_time : nest_out(i,3)+ref_time,2),'.');
end
title(strcat('befor & after time = ',num2str(ref_time/15),' sec'));
axis tight;

subplot(1,2,2);
plot([700 700],[pos(1,1) pos(end,1)],'r--');hold on;
plot([600 600],[pos(1,1) pos(end,1)],'r--');
plot(pos(:,3),pos(:,1),'.','color',([0.5 0.5 0.5]));
for i=1:size(nest_out,1)
    plot(pos(nest_out(i,3)-ref_time : nest_out(i,3)+ref_time,3),pos(nest_out(i,3)-ref_time : nest_out(i,3)+ref_time,1),'.');
end
title(strcat('befor & after time = ',num2str(ref_time/15),' sec'));
axis tight;

vt_out_name = strcat(vt_name(1:end-4),'_turn','.xls');
vt_out_fig_name = strcat(vt_name(1:end-4),'_turn','.jpg');
saveas(gcf,vt_out_fig_name);
writematrix(nest_out,vt_out_name);
