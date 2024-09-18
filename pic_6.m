
clear;clc;

filepath = ['E:\data\Wind_cmems\multi_sensor_data\'];
dirinput = dir(fullfile(filepath,'*.mat'));
playname = {dirinput.name};
for i=1:length(playname)
    year(i) = str2num(playname{i}(end-7:end-4));
end
[~ ,loc] = sort(year);
for i=1:length(loc)
    load([filepath,playname{loc(i)}]);
    if i==1
        u = wind_u;v = wind_v;
    else
        u = cat(3,u,wind_u);v = cat(3,v,wind_v);
    end
end
clear wind_u wind_v
for i=1997:2023
    for j=1:12
    date_uv((i-1997)*12+j) = datenum([sprintf('%d',i),sprintf('%02d',j)],'yyyymm');
    end
end
date_uv=date_uv(1:end-7);
lon = ncread('E:\data\Wind_cmems\ERS-2_SCAT_25_mm_1997_2000\GLO-WIND_L3-OBS_ERS-2_SCAT_25_ASC_19970101.nc','lon');
lat = ncread('E:\data\Wind_cmems\ERS-2_SCAT_25_mm_1997_2000\GLO-WIND_L3-OBS_ERS-2_SCAT_25_ASC_19970101.nc','lat');
lon = double(lon); lat = double(lat);
uv = sqrt(u.^2+v.^2);

filepath = ['E:\data\Wind_cmems\multi_sensor_data_wind_stress\wind_stress\'];
dirinput = dir(fullfile(filepath,'*.mat'));
playname = {dirinput.name};
for i=1:length(playname)
    year_s(i) = str2num(playname{i}(end-7:end-4));
end
[~ ,loc_s] = sort(year_s);
for i=1:length(loc_s)
    load([filepath,playname{loc_s(i)}]);
    if i==1
        u_str = u_stress;v_str = v_stress;
    else
        u_str = cat(3,u_str,u_stress);v_str = cat(3,v_str,v_stress);
    end
end
clear u_stress v_stress
for i=1997:2023
    for j=1:12
    date_str((i-1997)*12+j) = datenum([sprintf('%d',i),sprintf('%02d',j)],'yyyymm');
    end
end
date_str=date_str(1:end-7);


time = date_str(68:307);
ruo_a = 1.25;
ruo_w = 1023;
omiga = 7.29e-5;

tao_x_n = squeeze(nanmean(u_str(681:1040,379:382,68:307),2));
tao_x_s = squeeze(nanmean(u_str(681:1040,339:342,68:307),2));

et_n = nanmean(tao_x_n)./(2*omiga*ruo_w*sind(5));
et_s = nanmean(tao_x_s)./(2*omiga*ruo_w*sind(-5));
et_n = -1 * et_n;

et_all = abs(et_n) + abs(et_s);

[imf_2, residual_2] = emd(et_n,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(et_n,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(et_n,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(et_n,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(et_n,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:240
    et_n_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    et_n_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_n(:,1) = residual_2;
Tt_n(:,2) = residual_3;
Tt_n(:,3) = residual_4;
Tt_n(:,4) = residual_5;
Tt_n(:,5) = residual_6;

[imf_2, residual_2] = emd(et_s,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(et_s,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(et_s,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(et_s,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(et_s,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:240
    et_s_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    et_s_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_s(:,1) = residual_2;
Tt_s(:,2) = residual_3;
Tt_s(:,3) = residual_4;
Tt_s(:,4) = residual_5;
Tt_s(:,5) = residual_6;

[imf_2, residual_2] = emd(et_all,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(et_all,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(et_all,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(et_all,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(et_all,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:240
    et_all_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    et_all_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_a(:,1) = residual_2;
Tt_a(:,2) = residual_3;
Tt_a(:,3) = residual_4;
Tt_a(:,4) = residual_5;
Tt_a(:,5) = residual_6;

% Tt_n_yr = (et_n_m(end) - et_n_m(1))/242*12
Tt_n_m_2 = et_n_m - detrend(et_n_m);
Tt_n_yr_2 = (Tt_n_m_2(end) - Tt_n_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_n(:,i) - detrend(Tt_n(:,i));
    Tt_n_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_n = std(Tt_n_all);
% Tt_s_yr = (et_s_m(end) - et_s_m(1))/242*12
Tt_s_m_2 = et_s_m - detrend(et_s_m);
Tt_s_yr_2 = (Tt_s_m_2(end) - Tt_s_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_s(:,i) - detrend(Tt_s(:,i));
    Tt_s_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_s = std(Tt_s_all);
% Tt_w_yr = (et_all_m(end) - et_all_m(1))/242*12
Tt_w_m_2 = et_all_m - detrend(et_all_m);
Tt_w_yr_2 = (Tt_w_m_2(end) - Tt_w_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_a(:,i) - detrend(Tt_a(:,i));
    Tt_a_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_a = std(Tt_a_all);


%%
figure(11)
set(gcf,'pos',[2650 200 750 680])
set(gcf,'color',[1 1 1])

s1 = subplot(3,1,1)
hold on
f1 = fill([time';flipud(time')],[et_s_m*-10+2*et_s_std*-10;flipud(et_s_m*-10-2*et_s_std*-10)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time,et_s*-10,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time,et_s_m*-10,'color',[233 0 45]/256,'linewidth',1.2);

B=datestr(time(6+12:48:end),'yyyy');
set(gca, 'Box', 'off','layer','top','color','none', ...          % 边框
        'YAxisLocation','left',...
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-100 0],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-100:50:0,'yticklabel',-100:50:0);
ylabel({'Southward Ekman transport';'[Sv]'},'FontName','Arial',...
    'fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+285 0.85*100-100])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'pos',[0.12 0.702 0.77 0.258])
% t3 = text(time(150),-80,[sprintf('%5.3f',Tt_s_yr),' Sv/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time(120),-80,[sprintf('%5.2f',Tt_s_yr_2*240/12*10),'(±',sprintf('%4.2f',Tt_std_n*240/12*10),') Sv, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s2 = subplot(3,1,2)
hold on
f1 = fill([time';flipud(time')],[et_n_m*10+2*et_n_std*10;flipud(et_n_m*10-2*et_n_std*10)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time,et_n*10,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time,et_n_m*10,'color',[233 0 45]/256,'linewidth',1.2);

set(gca, 'Box', 'off','layer','top','color','none', ...         % 边框
         'YAxisLocation','right','linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[0 100],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',0:50:100,'yticklabel',0:50:100);
ylabel({'Northward Ekman transport';'[Sv]'},'Rotation',90,...
'FontName','Arial','fontsize',12,'fontweight','normal','position',[time(end)+430 50])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(end)-285 0.85*100])
set(gca,'pos',[0.12 0.394 0.77 0.258])
% t1 = text(time(174),85,[sprintf('%5.3f',Tt_n_yr),' Sv/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t2 = text(time(120),85,[sprintf('%5.2f',Tt_n_yr_2*240/12*10),'(±',sprintf('%4.2f',Tt_std_s*240/12*10),') Sv, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s3 = subplot(3,1,3)
hold on
f1 = fill([time';flipud(time')],[et_all_m*10+2*et_all_std*10;flipud(et_all_m*10-2*et_all_std*10)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time,et_all*10,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time,et_all_m*10,'color',[233 0 45]/256,'linewidth',1.2);
set(gca, 'Box', 'off','layer','top','color','none', ...              % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[0 150],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',0:75:150,'yticklabel',0:75:150);
ylabel({'Ekman divergence';'[Sv]'},'FontName','Arial','fontsize',12,'fontweight','normal')
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+285 0.85*150])
set(gca,'pos',[0.12 0.086 0.77 0.258])
% t5 = text(time(150-36),35,[sprintf('%5.3f',Tt_w_yr),' Sv/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t6 = text(time(120),120,[sprintf('%5.2f',Tt_w_yr_2*240/12*10),'(±',sprintf('%4.2f',Tt_std_a*240/12*10),') Sv, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

print(figure(11),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_6'],'-dpng','-r1200')

%%

n_t_yr = (et_n_m(end) - et_n_m(1))*10
s_t_yr = (et_s_m(end) - et_s_m(1))*10
d_t_yr = (et_all_m(end) - et_all_m(1))*10


