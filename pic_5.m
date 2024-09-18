
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

% 5°S-5°N(341:380);150°E-150°W(601:840)

u_m = squeeze(nanmean(nanmean(u(601:840,341:380,68:307))));

time = date_uv(68:307);

[imf_2, residual_2] = emd(u_m,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(u_m,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(u_m,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(u_m,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(u_m,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:240
    u_t_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    u_t_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end

Tt_u(:,1) = residual_2;
Tt_u(:,2) = residual_3;
Tt_u(:,3) = residual_4;
Tt_u(:,4) = residual_5;
Tt_u(:,5) = residual_6;

% aa = (u_t_m(end) - u_t_m(1))
% bb = (u_str_t_m(end) - u_str_t_m(1))/242*12;


file_uv = ['E:\data\OSCAR\data_mm_1_4_global\'];
dir_uv = dir(fullfile(file_uv,'*.mat'));
playname_uv = {dir_uv.name};
for i=1:length(playname_uv)
    year_uv(i,:) = playname_uv{i}(19:22);
end

for i=1:length(playname_uv)
    load([file_uv,playname_uv{i}]);
        u_c = u(341:380,681:1040,:);
    if i==1
        u_sec=u_c;
    else
        u_sec = cat(3,u_sec,u_c);
    end
end
clear u_c
lon_uv = ncread('E:\data\OSCAR\oscar_currents_final_19930101.nc','lon');
lat_uv = ncread('E:\data\OSCAR\oscar_currents_final_19930101.nc','lat');
for i=1:length(playname_uv)
    for j = 1:12
        time_uv((i-1)*12+j) = datenum([year_uv(i,:),sprintf('%02d',j),'01'],'yyyymmdd');
    end
end

find(time_uv==datenum('2002-08-01'))
find(time_uv==datenum('2022-09-01'))

u_c = squeeze(nanmean(nanmean(u_sec(:,:,116:355))));
time = time_uv(116:355);

[imf_2, residual_2] = emd(u_c,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(u_c,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(u_c,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(u_c,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(u_c,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:240
    u_c_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    u_c_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_c(:,1) = residual_2;
Tt_c(:,2) = residual_3;
Tt_c(:,3) = residual_4;
Tt_c(:,4) = residual_5;
Tt_c(:,5) = residual_6;

Tt_t_yr = (u_t_m(end) - u_t_m(1))/240*12
Tt_t_m_2 = u_t_m - detrend(u_t_m);
Tt_t_yr_2 = (Tt_t_m_2(end) - Tt_t_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_u(:,i) - detrend(Tt_u(:,i));
    Tt_u_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_u = std(Tt_u_all);

Tt_c_yr = (u_c_m(end) - u_c_m(1))/240*12
Tt_c_m_2 = u_c_m - detrend(u_c_m);
Tt_c_yr_2 = (Tt_c_m_2(end) - Tt_c_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_c(:,i) - detrend(Tt_c(:,i));
    Tt_c_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_c = std(Tt_c_all);



%%

figure(11)

set(gcf,'pos',[2650 200 770 500])
set(gcf,'color',[1 1 1])

s1=subplot(2,1,1)
hold on
f1 = fill([time';flipud(time')],[u_t_m+2*u_t_std;flipud(u_t_m-2*u_t_std)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time,u_m,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time,u_t_m,'color',[233 0 45]/256,'linewidth',1.2);

B=datestr(time(6+12:48:end),'yyyy');
set(gca, 'Box', 'off','layer','top','color','none', ...             % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-8 2],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-8:5:2,'yticklabel',-8:5:2);
ylabel({'Wind speed';'[m/s]'},'FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*10-8])
set(s1,'pos',[0.13 0.572 0.81 0.388])
% t1 = text(time(174),0.5,[sprintf('%5.3f',Tt_t_yr),' (m/s)/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t2 = text(time(120),-6.75,[sprintf('%5.2f',Tt_t_yr_2*240/12),'(±',sprintf('%4.2f',Tt_std_u*240/12),') m/s, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度


s2=subplot(2,1,2)
hold on
f1 = fill([time';flipud(time')],[u_c_m+2*u_c_std;flipud(u_c_m-2*u_c_std)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time,u_c,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time,u_c_m,'color',[233 0 45]/256,'linewidth',1.2);
set(gca, 'Box', 'off','layer','top','color','none', ...             % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-0.8 0.4],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-0.8:0.6:0.4,'yticklabel',-0.8:0.6:0.4);
ylabel({'Zonal SEC';'[m/s]'},'FontName','Arial','fontsize',12,'fontweight','normal')
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*1.2-0.8])
set(s2,'pos',[0.13 0.11 0.81 0.388])
% t3 = text(time(174),0.3,[sprintf('%6.4f',Tt_c_yr),' (m/s)/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time(120),-0.65,[sprintf('%6.2f',Tt_c_yr_2*240/12),'(±',sprintf('%4.2f',Tt_std_c*240/12),') m/s, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

print(figure(11),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_5'],'-dpng','-r1200')

%%

u_t_yr = (u_t_m(end) - u_t_m(1))
u_str_yr = (u_str_t_m(end) - u_str_t_m(1))
u_c_yr = (u_c_m(end) - u_c_m(1))






