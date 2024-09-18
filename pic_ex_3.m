clear;clc;

filepath_chl=['E:\data\chl_9km_monthly\'];

dirOutput = dir(fullfile(filepath_chl,'*.nc'));
playname={dirOutput.name};
for i=1:240;
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

load('E:\matlab_code\Tropical_Pacific\tropical_high_chl_region2.mat');
highc_re=[line_west;flipud(line_east);line_west(1,:)];


R=6371;
tro_pac_in=inpolygon(X,Y,highc_re(:,1),highc_re(:,2));

S_all = sum(sum((2*pi*R/360*(lat(1)-lat(2)))*(2*pi*R/360*(lon(2)-lon(1))*cosd(Y)).*tro_pac_in));
S_sea = sum(sum((2*pi*R/360*(lat(1)-lat(2)))*(2*pi*R/360*(lon(2)-lon(1))*cosd(Y)).*tro_pac_in));
s_fac=(2*pi*R/360*(lat(1)-lat(2)))*(2*pi*R/360*(lon(2)-lon(1))*cosd(Y));
ratio_sea_all=S_sea/S_all;

for i=1:240
    in_chl10=zeros(2400,1200);
    chl_temp=chl(:,:,i);
    for j=1:2400
        for k=1:1200
            if chl_temp(j,k)>0.1;
                in_chl10(j,k)=1;
            end
        end
    end
    S_highchl10_m(i)=sum(sum(s_fac.*tro_pac_in.*in_chl10'));
end


[imf_2, residual_2] = emd(S_highchl10_m,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(S_highchl10_m,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(S_highchl10_m,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(S_highchl10_m,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(S_highchl10_m,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:240
    residual_std_1(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    residual_m_1(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_1(:,1) = residual_2;
Tt_1(:,2) = residual_3;
Tt_1(:,3) = residual_4;
Tt_1(:,4) = residual_5;
Tt_1(:,5) = residual_6;

filename = ['E:\data\OCCCI\fv6_4km_1m\ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-199709-fv6.0.nc'];
lon = ncread(filename,'lon');
lat = ncread(filename,'lat');
lon = cat(1,lon(4321:end),lon(1:4320)+360);
[X,Y] = meshgrid(lon(2401:7200),lat(961:3360));

input_path = ['E:\data\OCCCI\fv6_4km_1m\data_mm_4km_global\'];
dirinput = dir(fullfile(input_path,'*.mat'));
playname = {dirinput.name};

for i=1:length(playname)
    date_time(i,:) = playname{i}(18:23);
end
time_c = datenum(date_time,'yyyymm');

load('E:\matlab_code\Tropical_Pacific\tropical_high_chl_region2.mat');
highc_re = [line_west;flipud(line_east);line_west(1,:)];
R = 6371;
tro_pac_in = inpolygon(X,Y,highc_re(:,1),highc_re(:,2));
S_all = sum(sum((2*pi*R/360*(lat(1)-lat(2)))*(2*pi*R/360*(lon(2)-lon(1))*cosd(Y)).*tro_pac_in));
s_fac = (2*pi*R/360*(lat(1)-lat(2)))*(2*pi*R/360*(lon(2)-lon(1))*cosd(Y));

for i=1:length(playname)
    load([input_path,playname{i}]);
    chl_temp = cat(1,chlor_a(4321:end,:),chlor_a(1:4320,:));
    chl_temp = chl_temp(2401:7200,961:3360);
    
    in_chl10 = zeros(4800,2400);
    for j = 1:4800
        for k = 1:2400
            if chl_temp(j,k)>0.1;
                in_chl10(j,k) = 1;
            end
        end
    end
    S_highchl10(i) = sum(sum(s_fac.*tro_pac_in.*in_chl10'));
end

time_c_2 = time_c(1:312);
S_hct = S_highchl10(1:312);
%%
[imf, residual] = emd(S_hct);
[imf_2, residual_2] = emd(S_hct,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(S_hct,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(S_hct,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(S_hct,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(S_hct,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:312
    residual_std_2(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:312
    residual_m_2(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_2(:,1) = residual_2;
Tt_2(:,2) = residual_3;
Tt_2(:,3) = residual_4;
Tt_2(:,4) = residual_5;
Tt_2(:,5) = residual_6;

time_c_3 = time_c(37:312);
S_hct_2 = S_highchl10(37:312);

% [imf, residual] = emd(S_hct);
[imf_2, residual_2] = emd(S_hct_2,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(S_hct_2,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(S_hct_2,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(S_hct_2,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(S_hct_2,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

for i=1:276
    residual_std_3(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
    residual_m_3(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
Tt_3(:,1) = residual_2;
Tt_3(:,2) = residual_3;
Tt_3(:,3) = residual_4;
Tt_3(:,4) = residual_5;
Tt_3(:,5) = residual_6;

Tt_1_yr = (residual_m_1(end) - residual_m_1(1))/240*12
Tt_1_m_2 = residual_m_1 - detrend(residual_m_1);
Tt_1_yr_2 = (Tt_1_m_2(end) - Tt_1_m_2(1))/240*12
for i=1:5
    Tt_m_s = Tt_1(:,i) - detrend(Tt_1(:,i));
    Tt_1_all(i) = (Tt_m_s(end) - Tt_m_s(1))/240*12;
end
Tt_1_std = std(Tt_1_all);

Tt_2_yr = (residual_m_2(end) - residual_m_2(1))/312*12
Tt_2_m_2 = residual_m_2 - detrend(residual_m_2);
Tt_2_yr_2 = (Tt_2_m_2(end) - Tt_2_m_2(1))/312*12
for i=1:5
    Tt_m_s = Tt_2(:,i) - detrend(Tt_2(:,i));
    Tt_2_all(i) = (Tt_m_s(end) - Tt_m_s(1))/312*12;
end
Tt_2_std = std(Tt_2_all)

Tt_3_yr = (residual_m_3(end) - residual_m_3(1))/276*12
Tt_3_m_2 = residual_m_3 - detrend(residual_m_3);
Tt_3_yr_2 = (Tt_3_m_2(end) - Tt_3_m_2(1))/276*12
for i=1:5
    Tt_m_s = Tt_3(:,i) - detrend(Tt_3(:,i));
    Tt_3_all(i) = (Tt_m_s(end) - Tt_m_s(1))/276*12;
end
Tt_3_std = std(Tt_3_all)

%%

figure(11)
set(gcf,'pos',[250 200 750 680])
set(gcf,'color',[1 1 1])

s1 = subplot(3,1,1)
hold on
f1 = fill([time_chl';flipud(time_chl')],[residual_m_1+2*residual_std_1;flipud(residual_m_1-2*residual_std_1)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time_chl,S_highchl10_m,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time_chl,residual_m_1,'color',[233 0 45]/256,'linewidth',1.2)

B=datestr(time_c(29:48:end),'yyyy');
set(gca,'Xtick',time_c(29:48:end),'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',0:2e7:4e7,'yticklabel',0:2:4,'fontsize',12,'fontweight','normal');
ylabel({'CRT area from MODIS';'[10^7 km^2]'},'fontsize',12,'fontweight','normal')
xlim([time_c(1) time_c(end)])
ylim([0 4e7])
% title('Time series of HCT area in the Tropical Pacific','fontsize',15,'fontweight','bold');
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1]);           % 坐标轴颜色
set(gca,'layer','top')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_c(1)+340 0.85*4e7])
set(gca,'pos',[0.12 0.7 0.77 0.26])
% t1 = text(time_chl(150-48),1.3e7,[sprintf('%4.2f',8.09),'×10^4 km^2/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t2 = text(time_chl(50),1e7,['1.87(±0.82) × 10^5 km^2/yr, {\itP} < 0.05'],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s2 = subplot(3,1,2)
hold on
f1 = fill([time_c;flipud(time_c)],[residual_m_2+residual_std_2;flipud(residual_m_2-residual_std_2)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time_c,S_highchl10,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time_c,residual_m_2,'color',[233 0 45]/256,'linewidth',1.2)
xlim([time_c(1) time_c(end)])
ylim([0 4e7])
set(gca,'Xtick',time_c(29:48:end),'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',0:2e7:4e7,'yticklabel',0:2:4,'fontsize',12,'fontweight','normal');
ylabel({'CRT area from OC-CCI';'[10^7 km^2]'},'Rotation',90,...
    'fontsize',12,'fontweight','normal','position',[time_c(end)+400 2e7])
% title('Time series of HCT area in the Tropical Pacific','fontsize',15,'fontweight','bold');
set(gca, 'Box', 'off','color','none', ...            % 边框
        'YAxisLocation','right',...
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1]);           % 坐标轴颜色
set(gca,'layer','top')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_c(end)-370 0.85*4e7])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'pos',[0.12 0.3925 0.77 0.26])
% t3 = text(time_chl(150-48),1.3e7,[sprintf('%4.2f',4.76),'×10^5 km^2/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time_chl(50),1.6e7,['4.57(±0.72) × 10^5 km^2/yr, {\itP} < 0.05'],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s3 = subplot(3,1,3)
hold on
f1 = fill([time_c_3;flipud(time_c_3)],[residual_m_3+residual_std_3;flipud(residual_m_3-residual_std_3)],[.7 .7 .7])
set(f1,'EdgeColor','none','FaceAlpha',0.5)
l1 = plot(time_c_3,S_hct_2,'color',[2 38 62]/256,'linewidth',1.3)
l2 = plot(time_c_3,residual_m_3,'color',[233 0 45]/256,'linewidth',1.2)

xlim([time_c(1) time_c(end)])
ylim([0 4e7])
set(gca,'Xtick',time_c(29:48:end),'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',0:2e7:4e7,'yticklabel',0:2:4,'fontsize',12,'fontweight','normal');
ylabel({'CRT area from OC-CCI';'[10^7 km^2]'},'fontsize',12,'fontweight','normal')
% title('Time series of HCT area in the Tropical Pacific','fontsize',15,'fontweight','bold');
set(gca, 'Box', 'off','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1]);           % 坐标轴颜色
set(gca,'layer','top')
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_c(1)+340 0.85*4e7])
set(gca,'pos',[0.12 0.085 0.77 0.26])
% t5 = text(time_chl(150-48),1.3e7,[sprintf('%4.2f',4.27),'×10^4 km^2/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t6 = text(time_chl(50),1.6e7,['7.58(±1.82) × 10^4 km^2/yr, {\itP} < 0.05'],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

print(figure(11),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_3_2'],'-dpng','-r1200') 

%%
t1 = (residual_m_1(end)-residual_m_1(1))/length(residual_m_1)*12
t2 = (residual_m_2(end)-residual_m_2(1))/length(residual_m_2)*12
t3 = (residual_m_3(end)-residual_m_3(1))/length(residual_m_3)*12



