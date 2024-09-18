
clear;clc;

filepath_chl=['E:\data\chl_9km_monthly\'];

dirOutput = dir(fullfile(filepath_chl,'*.nc'));
playname={dirOutput.name};
for i=1:242;
    dirname=playname{i};
    chl_temp=ncread([filepath_chl,dirname],'chlor_a');
    chl_temp=cat(1,chl_temp(2161:end,:),chl_temp(1:2160,:));
    chl_temp=chl_temp(1201:3600,481:1680);
    chl(:,:,i)=chl_temp;
    time_chl(i)=datenum(dirname(12:19),'yyyymmdd');
end
%lon:100E-300E:1201:3600
lon=ncread([filepath_chl,playname{1}],'lon');
lon=cat(1,lon(2161:end),lon(1:2160)+360);
lon=lon(1201:3600);
%lat:50S-50N:481:1680
lat=ncread([filepath_chl,playname{1}],'lat');
lat=lat(481:1680);
[X,Y]=meshgrid(lon,lat);
time_chl = time_chl(1:240);
%%

chl_lat = squeeze(nanmean(chl(1441:1560,361:840,:),1));
lat_loc = flipud(lat(361:600));

for i=1:240
    chl_temp = chl_lat(1:240,i);
    chl_temp = flipud(chl_temp)-0.1;
    lon_temp = find(chl_temp<0);
    if length(lon_temp)<13
        lat_n(i) = nan;
    else
        lon_temp2 = lon_temp(13:end) - lon_temp(1:end-12);
        lon_temp3 = find(lon_temp2==12);
        if length(lon_temp3)>0
            lat_n(i) = lat_loc(lon_temp(lon_temp3(1)));
        else
            lat_n(i) = nan;
        end
    end
end

time_interp = time_chl(~isnan(lat_n));
lat_interp = lat_n(~isnan(lat_n));
lat_n_interp = interp1(time_interp,lat_interp,time_chl);

[imf_2, residual_2] = emd(lat_n_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lat_n_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lat_n_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lat_n_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lat_n_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

St_n(:,1) = imf_2(:,1) + imf_2(:,2);
St_n(:,2) = imf_3(:,1) + imf_3(:,2);
St_n(:,3) = imf_4(:,1) + imf_4(:,2);
St_n(:,4) = imf_5(:,1) + imf_5(:,2);
St_n(:,5) = imf_6(:,1) + imf_6(:,2);

It_n(:,1) = imf_2(:,3) + imf_2(:,4) + imf_2(:,5) + imf_2(:,6);
It_n(:,2) = imf_3(:,3) + imf_3(:,4) + imf_3(:,5) + imf_3(:,6);
It_n(:,3) = imf_4(:,3) + imf_4(:,4) + imf_4(:,5) + imf_4(:,6);
It_n(:,4) = imf_5(:,3) + imf_5(:,4) + imf_5(:,5) + imf_5(:,6);
It_n(:,5) = imf_6(:,3) + imf_6(:,4) + imf_6(:,5) + imf_6(:,6);

Tt_n(:,1) = residual_2;
Tt_n(:,2) = residual_3;
Tt_n(:,3) = residual_4;
Tt_n(:,4) = residual_5;
Tt_n(:,5) = residual_6;

for i=1:240
    St_n_m(i) = nanmean(St_n(i,:));
    It_n_m(i) = nanmean(It_n(i,:));
    Tt_n_m(i) = nanmean(Tt_n(i,:));
    St_n_std(i) = std(St_n(i,:));
    It_n_std(i) = std(It_n(i,:));
    Tt_n_std(i) = std(Tt_n(i,:));
end

lat_loc2 = lat(601:840);

for i=1:240
    chl_temp = chl_lat(241:480,i);
    chl_temp = chl_temp-0.1;
    lon_temp = find(chl_temp<0);
    if length(lon_temp)<13
        lat_s(i) = nan;
    else
        lon_temp2 = lon_temp(13:end) - lon_temp(1:end-12);
        lon_temp3 = find(lon_temp2==12);
        if length(lon_temp3)>0
            lat_s(i) = lat_loc2(lon_temp(lon_temp3(1)));
        else
            lat_s(i) = nan;
        end
    end
end

time_interp = time_chl(~isnan(lat_s));
lat_interp = lat_s(~isnan(lat_s));
lat_s_interp = interp1(time_interp,lat_interp,time_chl);
lat_s_interp(105) = lat_loc2(147);
lat_s_interp(115) = lat_loc2(162);
lat_s_interp(116) = lat_loc2(156);

[imf_2, residual_2] = emd(lat_s_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lat_s_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lat_s_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lat_s_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lat_s_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

St_s(:,1) = imf_2(:,1) + imf_2(:,2);
St_s(:,2) = imf_3(:,1) + imf_3(:,2);
St_s(:,3) = imf_4(:,1) + imf_4(:,2);
St_s(:,4) = imf_5(:,1) + imf_5(:,2);
St_s(:,5) = imf_6(:,1) + imf_6(:,2);

It_s(:,1) = imf_2(:,3) + imf_2(:,4);
It_s(:,2) = imf_3(:,3) + imf_3(:,4) + imf_3(:,5);
It_s(:,3) = imf_4(:,3) + imf_4(:,4) + imf_4(:,5);
It_s(:,4) = imf_5(:,3) + imf_5(:,4) + imf_5(:,5);
It_s(:,5) = imf_6(:,3) + imf_6(:,4) + imf_6(:,5);

Tt_s(:,1) = residual_2;
Tt_s(:,2) = residual_3;
Tt_s(:,3) = residual_4;
Tt_s(:,4) = residual_5;
Tt_s(:,5) = residual_6;

for i=1:240
    St_s_m(i) = nanmean(St_s(i,:));
    It_s_m(i) = nanmean(It_s(i,:));
    Tt_s_m(i) = nanmean(Tt_s(i,:));
    St_s_std(i) = std(St_s(i,:));
    It_s_std(i) = std(It_s(i,:));
    Tt_s_std(i) = std(Tt_s(i,:));
end

chl_equ_2 = squeeze(nanmean(chl(481:1440,589:612,:),2));
lon_loc = flipud(lon(481:1440));
for i=1:240
    chl_temp = chl_equ_2(:,i);
    chl_temp = flipud(chl_temp)-0.1;
    lon_temp = find(chl_temp<0);
    if length(lon_temp)<4
        lon_2(i) = nan;
    else
        lon_temp2 = lon_temp(3:end) - lon_temp(1:end-2);
        lon_temp3 = find(lon_temp2==2);
        if length(lon_temp3)>0
            lon_2(i) = lon_loc(lon_temp(lon_temp3(1)));
        else
            lon_2(i) = nan;
        end
    end
end
time_interp = time_chl(~isnan(lon_2));
lon_interp = lon_2(~isnan(lon_2));
lon_2_interp_1 = interp1(time_interp,lon_interp,time_chl,'linear');
lon_2_interp = lon_2_interp_1(1:237);

[imf_2, residual_2] = emd(lon_2_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lon_2_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lon_2_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lon_2_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lon_2_interp,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

St_w(:,1) = imf_2(:,1) + imf_2(:,2);
St_w(:,2) = imf_3(:,1) + imf_3(:,2);
St_w(:,3) = imf_4(:,1) + imf_4(:,2);
St_w(:,4) = imf_5(:,1) + imf_5(:,2);
St_w(:,5) = imf_6(:,1) + imf_6(:,2);

It_w(:,1) = imf_2(:,3) + imf_2(:,4) + imf_2(:,5) + imf_2(:,6);
It_w(:,2) = imf_3(:,3) + imf_3(:,4) + imf_3(:,5) + imf_3(:,6);
It_w(:,3) = imf_4(:,3) + imf_4(:,4) + imf_4(:,5) + imf_4(:,6);
It_w(:,4) = imf_5(:,3) + imf_5(:,4) + imf_5(:,5) + imf_5(:,6);
It_w(:,5) = imf_6(:,3) + imf_6(:,4) + imf_6(:,5) + imf_6(:,6);

Tt_w(:,1) = residual_2;
Tt_w(:,2) = residual_3;
Tt_w(:,3) = residual_4;
Tt_w(:,4) = residual_5;
Tt_w(:,5) = residual_6;

for i=1:237
    St_w_m(i) = nanmean(St_w(i,:));
    It_w_m(i) = nanmean(It_w(i,:));
    Tt_w_m(i) = nanmean(Tt_w(i,:));
    St_w_std(i) = std(St_w(i,:));
    It_w_std(i) = std(It_w(i,:));
    Tt_w_std(i) = std(Tt_w(i,:));
end


Tt_n_yr = (Tt_n_m(end) - Tt_n_m(1))/240*12
Tt_n_m_2 = Tt_n_m - detrend(Tt_n_m);
Tt_n_yr_2 = (Tt_n_m_2(end) - Tt_n_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_n(:,i) - detrend(Tt_n(:,i));
    Tt_n_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_n = std(Tt_n_all);

Tt_s_yr = (Tt_s_m(end) - Tt_s_m(1))/240*12
Tt_s_m_2 = Tt_s_m - detrend(Tt_s_m);
Tt_s_yr_2 = (Tt_s_m_2(end) - Tt_s_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_s(:,i) - detrend(Tt_s(:,i));
    Tt_s_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_s = std(Tt_s_all);

Tt_w_yr = (Tt_w_m(end) - Tt_w_m(1))/237*12
Tt_w_m_2 = Tt_w_m - detrend(Tt_w_m);
Tt_w_yr_2 = (Tt_w_m_2(end) - Tt_w_m_2(1))/237*12
for i=1:5
    Tt_dt = Tt_w(:,i) - detrend(Tt_w(:,i));
    Tt_w_all(i) = (Tt_dt(end) - Tt_dt(1))/237*12;
end
Tt_std_w = std(Tt_w_all);

%%

figure(11)
set(gcf,'pos',[2650 200 750 680])
set(gcf,'color',[1 1 1])

s1 = subplot(3,1,1)
hold on
f3 = fill([time_chl(1:237)';flipud(time_chl(1:237)')],[Tt_w_m'+2*Tt_w_std';flipud(Tt_w_m'-2*Tt_w_std')],[.7 .7 .7])
l3 = plot(time_chl(1:237),lon_2_interp,'color',[2 38 62]/256,'linewidth',1.3)
k3 = plot(time_chl(1:237),Tt_w_m,'color',[233 0 45]/256,'linewidth',1.2)
set(f3,'EdgeColor','none','FaceAlpha',0.5)
B=datestr(time_chl(6+12:48:end),'yyyy');
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[140 220],...
         'Xtick',time_chl(6+12:48:end),'xticklabel',B,...
         'Ytick',140:40:220,'yticklabel',['140°E';'180° ';'140°W']);
ylabel('Western boundary','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(1)+265 0.85*80+140])
set(gca,'pos',[0.12 0.705 0.77 0.26])
% t1 = text(time_chl(174),200,[sprintf('%4.2f',Tt_w_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t2 = text(time_chl(52),200,[sprintf('%4.2f',Tt_w_yr_2),'(±',sprintf('%4.2f',Tt_std_w),') °/yr, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s2 = subplot(3,1,2)
hold on
f1 = fill([time_chl';flipud(time_chl')],[Tt_n_m'+2*Tt_n_std';flipud(Tt_n_m'-2*Tt_n_std')],[.7 .7 .7])
l1 = plot(time_chl,lat_n_interp,'color',[2 38 62]/256,'linewidth',1.3)
k1 = plot(time_chl,Tt_n_m,'color',[233 0 45]/256,'linewidth',1.2)
set(f1,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
    'YAxisLocation','right',...
         'linewidth', 1.2,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[0 20],...
         'Xtick',time_chl(6+12:48:end),'xticklabel',B,...
         'Ytick',0:10:20,'yticklabel',[' 0° ';'10°N';'20°N']);
ylabel('Northern boundary','Rotation',90,'FontName','Arial',...
    'fontsize',12,'fontweight','normal','position',[time_chl(end)+620 10])
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-400 0.85*20])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'pos',[0.12 0.395 0.77 0.26])
% t3 = text(time_chl(150),17,[sprintf('%5.3f',Tt_n_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time_chl(46),4,[sprintf('%5.3f',Tt_n_yr_2),'(±',sprintf('%4.2f',Tt_std_n),') °/yr, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s3 = subplot(3,1,3)
hold on
f2 = fill([time_chl';flipud(time_chl')],[Tt_s_m'+2*Tt_s_std';flipud(Tt_s_m'-2*Tt_s_std')],[.7 .7 .7])
l2 = plot(time_chl,lat_s_interp,'color',[2 38 62]/256,'linewidth',1.3)
k2 = plot(time_chl,Tt_s_m,'color',[233 0 45]/256,'linewidth',1.2)
set(f2,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.2,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[-20 0],...
         'Xtick',time_chl(6+12:48:end),'xticklabel',B,...
         'Ytick',-20:10:0,'yticklabel',['20°S';'10°S';' 0° ']);
ylabel('Southern boundary','FontName','Arial','fontsize',12,'fontweight','normal')
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(1)+265 0.85*20-20])
set(gca,'pos',[0.12 0.085 0.77 0.26])
% t5 = text(time_chl(150-36),-5,[sprintf('%5.3f',Tt_s_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t6 = text(time_chl(46),-8.8,[sprintf('%5.3f',Tt_s_yr_2),'(±',sprintf('%4.2f',Tt_std_s),') °/yr, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

print(figure(11),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_3'],'-dpng','-r1200')

%%

Tt_n_yr = (Tt_n_m(end) - Tt_n_m(1))/242*12
Tt_s_yr = (Tt_s_m(end) - Tt_s_m(1))/242*12
Tt_w_yr = (Tt_w_m(end) - Tt_w_m(1))/237*12




