function TraceAnal_any

f = figure('units','normalized','outerposition',[0 0 1 1]);

p = uipanel(f,'Position',[0.75 0.1 0.2 0.05]);
c = uicontrol(p,'Style','pushbutton','Position',[220 10 100 30]);
c.String = 'Close';
c.Callback = @plotButtonPushed;

% p0 = uipanel(f,'Position',[0.8 0.1 0.05 0.05]);
c0 = uicontrol(p,'Style','pushbutton','Position',[50 10 100 30]);
c0.String = 'Open';
c0.Callback = @plotButtonPushed2;
file='';
pos = 0;
pos_sp_y=0;
pos_sp_x=0;
grayColor = [.7 .7 .7];
x0 = 0;


    function plotButtonPushed2(src,event)
        [file,path] = uigetfile('*.csv');
        vt_name = file;%'H18_010422_VT1.nvt';

        pos = cal_vt_any(vt_name);

        x0 = pos(1,1):10000:pos(end,1);        
        pos_sp_y = interp1(pos(:,1),pos(:,2),x0);
        pos_sp_x = interp1(pos(:,1),pos(:,3),x0);

        subplot(1,4,1);
        plot(pos_sp_x,x0,'r.');hold on;
        plot(pos(:,3),pos(:,1),'.','Color', grayColor);
        subplot(1,4,2);
        plot(x0,pos_sp_y,'r.');hold on;
        plot(pos(:,1),pos(:,2),'.','Color', grayColor);
        subplot(1,4,3);
        plot(pos_sp_x,pos_sp_y,'r.');hold on;
        plot(pos(:,3),pos(:,2),'.','Color', grayColor);
        subplot(1,4,3);
    end

p2 = uipanel(f,'Position',[0.75 0.15 0.2 0.05],'Title','Windows Size');
c1 = uicontrol(p2,'Style','slider','position',[10,10,350,20],'SliderStep',[1/500 .05]);
c1.Value = 1;
c1.Callback = @sliderMoving;

p3 = uipanel(f,'Position',[0.75 0.2 0.2 0.05],'Title','Time location');
c2 = uicontrol(p3,'Style','slider','position',[10,10,350,20]);
c2.Value = 0;
c2.Callback = @sliderMoving2;

time_visible = size(pos,1);



    function plotButtonPushed(src,event)
        close all;
    end

    function sliderMoving(event,cg)
        if(event.Value==0)
            event.Value=0.01;
        end
        len = size(pos,1);
        subplot(1,4,1);
        ylim([pos(1,1) pos(floor(len*event.Value),1)]);xlim([200 1000]);
        subplot(1,4,2);
        xlim([pos(1,1) pos(floor(len*event.Value),1)]);ylim([0 900]);
        subplot(1,4,3);hold off;
        plot(pos(1:(floor(len*event.Value)),3),pos(1:(floor(len*event.Value)),2),'.','Color', grayColor);ylim([0 900]);xlim([200 1000]);
        time_visible = (floor(len*event.Value));
        win_size = floor(size(pos,1)/time_visible);
        c2.Value = 0;
        c2.Max = win_size;
        c2.SliderStep = [1/win_size 10/win_size];
    end

    function sliderMoving2(event,cg)
        if(event.Value==0)
            event.Value=1;
        end
        event.Value;
%         pos(floor((event.Value-1)*time_visible),1)
%         pos(floor((event.Value)*time_visible),1)
        
        subplot(1,4,1);
        ylim([pos(1+floor((event.Value-1)*time_visible),1) pos(1+floor((event.Value)*time_visible),1)]);
        subplot(1,4,2);
        xlim([pos(1+floor((event.Value-1)*time_visible),1) pos(1+floor((event.Value)*time_visible),1)]);
        subplot(1,4,3);hold off;
        plot(pos(1+floor((event.Value-1)*time_visible):1+floor((event.Value)*time_visible),3),pos(1+floor((event.Value-1)*time_visible):1+floor((event.Value)*time_visible),2),'.','Color', grayColor);ylim([0 900]);xlim([200 1000]);
    end

p4 = uipanel(f,'Position',[0.75 0.3 0.2 0.07],'Title','Save (filename)');
c4 = uicontrol(p4,'Style','pushbutton','Position',[250 10 100 30]);
c4.String = 'Save';
c4.Callback = @plotButtonPushed3;
c5 = uicontrol(p4,'Style','edit','Position',[30 10 200 30]);
c5.String = 'Nest_out';

    function plotButtonPushed3(src,event)
        vt_out_name = strcat(file,'_',c5.String,'.xls')
        writematrix(list,vt_out_name);

        dcm.removeAllDataCursors;
        c8.String='';
        c9.String='';
        c10.String='';
    end

p5 = uipanel(f,'Position',[0.75 0.4 0.2 0.07],'Title','Choose timepoints');
c6 = uicontrol(p5,'Style','pushbutton','Position',[50 10 100 30],'String','Rec on');
c6.Callback = @plotButtonPushed4;
c7 = uicontrol(p5,'Style','pushbutton','Position',[250 10 100 30],'String','Rec off');
c7.Callback = @plotButtonPushed5;

dcm = datacursormode(f);
dcm.SnapToDataVertex = 'on';
vals = null(1,1);
list = null(1,1);

    function plotButtonPushed4(src,event)
        dcm.Enable = 'on';
    end

    function plotButtonPushed5(src,event)
        vals = getCursorInfo(dcm);

        for i=1:size(vals,2)
            temp(i,:) = [vals(i).Position vals(i).DataIndex];
        end

        list = sortrows(temp,3,'ascend');

        c8.String = mat2str(list(:,1));
        c9.String = mat2str(list(:,2));
        c10.String = mat2str(list(:,3));
    end

p6 = uipanel(f,'Position',[0.75 0.5 0.2 0.3],'Title','List');
c8 = uicontrol(p6,'Style','text','Position',[10 10 30 350],'String','test');
c9 = uicontrol(p6,'Style','text','Position',[100 10 30 350],'String','test');
c10 = uicontrol(p6,'Style','text','Position',[200 10 30 350],'String','test');

p7 = uipanel(f,'Position',[0.75 0.82 0.2 0.1],'Title','Guide line');
c11 = uicontrol(p7,'Style','edit','Position',[10 10 70 30],'String','310');
c12= uicontrol(p7,'Style','edit','Position',[100 10 70 30],'String','235');
c13 = uicontrol(p7,'Style','pushbutton','Position',[200 10 100 30],'String','Plot');
c13.Callback = @plotButtonPushed6;

    function plotButtonPushed6(src,event)
        guide_x=str2double(c11.String);
        guide_y=str2double(c12.String);
        subplot(1,4,1);hold on;
%         plot([guide_x guide_x],[pos(1,1) pos(end,1)],'Color',[1 0.6 0.6]);
        xline(guide_x,'--')
        subplot(1,4,2);hold on;
%         plot([pos(1,1) pos(end,1)],[guide_y guide_y],'Color',[1 0.6 0.6]);
        yline(guide_y,'--')
        subplot(1,4,3);hold on;
        xline(guide_x,'--')
        yline(guide_y,'--')
    end

c14 = uicontrol(p7,'Style','text','Position',[10 40 30 30],'String','X');
c14 = uicontrol(p7,'Style','text','Position',[100 40 30 30],'String','Y');


end
