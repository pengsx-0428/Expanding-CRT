clear;clc;

load('E:\进展\文章相关\EEP_HCT\nc_reviewer\code\S_crt.mat');
S_crt = S_crt(1:240);
time = time(1:240);
% S_ext = Extend_sig_v2(S_crt,{'asymw'},605);

colortable=textread('D:\Matlab\anzhuang\bin\colorbar_NCL\MPL_rainbow.txt');
colortable_l = colortable(round(2:126/21:127),:);


%%

[imf_2, residual_2] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

St(:,1) = imf_2(:,1) + imf_2(:,2);
St(:,2) = imf_3(:,1) + imf_3(:,2);
St(:,3) = imf_4(:,1) + imf_4(:,2);
St(:,4) = imf_5(:,1) + imf_5(:,2);
St(:,5) = imf_6(:,1) + imf_6(:,2);

It(:,1) = imf_2(:,3) + imf_2(:,4) + imf_2(:,5) + imf_2(:,6);
It(:,2) = imf_3(:,3) + imf_3(:,4) + imf_3(:,5);
It(:,3) = imf_4(:,3) + imf_4(:,4) + imf_4(:,5);
It(:,4) = imf_5(:,3) + imf_5(:,4) + imf_5(:,5);
It(:,5) = imf_6(:,3) + imf_6(:,4) + imf_6(:,5) + imf_6(:,6);

% Et(:,1) = imf_2(:,3) + imf_2(:,4) + imf_2(:,5);
% Et(:,2) = imf_3(:,3) + imf_3(:,4) + imf_3(:,5);
% Et(:,3) = imf_4(:,3) + imf_4(:,4) + imf_4(:,5);
% Et(:,4) = imf_5(:,3) + imf_5(:,4) + imf_5(:,5);
% Et(:,5) = imf_6(:,3) + imf_6(:,4) + imf_6(:,5);

Tt(:,1) = residual_2;
Tt(:,2) = residual_3;
Tt(:,3) = residual_4;
Tt(:,4) = residual_5;
Tt(:,5) = residual_6;

for i=1:240
    St_m(i) = nanmean(St(i,:));
    It_m(i) = nanmean(It(i,:));
    Tt_m(i) = nanmean(Tt(i,:));
    St_std(i) = std(St(i,:));
    It_std(i) = std(It(i,:));
    Tt_std(i) = std(Tt(i,:));
end
% for i=1:240
%     Et_m(i) = nanmean(Et(i,:));
%     Et_std(i) = std(Et(i,:));
% end
S_unit = S_crt - It_m';
% S_unet = S_crt - Et_m';
for i=1:12
    S_first(i) = mean(S_unit(i:12:48+i));
    S_end(i) = mean(S_unit(180+i:12:228+i));
end
S_first = [S_first(6:12),S_first(1:5)];
S_end = [S_end(6:12),S_end(1:5)];

%%
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

load('E:\data\chl_9km_monthly\occci_in.mat');
in_temp=cat(2,in_temp(:,2161:end),in_temp(:,1:2160));
sea_in=~in_temp(481:1680,1201:3600);
x=x(481:1680,1201:3600);y=y(481:1680,1201:3600);
load('E:\matlab_code\Tropical_Pacific\tropical_high_chl_region2.mat');
line_west2 = [131.7,20;135.6,10.1;151,0;160,-4;180,-8;205,-20];
highc_re=[line_west2;flipud(line_east);line_west2(1,:)];
tro_pac_in=inpolygon(X,Y,highc_re(:,1),highc_re(:,2));

%%
month_1 = [1:60];
month_2 = [181:240];
chl_1 = nanmean(chl(:,:,month_1),3);
chl_2 = nanmean(chl(:,:,month_2),3);

chl_log_36 = zeros(2400,1200);
chl_log_215 = zeros(2400,1200);

for i = 1:2400
    for j = 1:1200
        if chl_1(i,j)>=0.1
            chl_log_36(i,j) = 1;
        end
        if chl_2(i,j)>=0.1
            chl_log_215(i,j) = 1;
        end
    end
end

chl_first = zeros(2400,1200);
chl_both = zeros(2400,1200);
chl_end = zeros(2400,1200);

for i=1:2400
    for j=1:1200
        if (chl_log_36(i,j)==1 & chl_log_215(i,j)==0)
            chl_first(i,j) = 1;
        elseif (chl_log_36(i,j)==1 & chl_log_215(i,j)==1)
            chl_both(i,j) = 1;
        elseif (chl_log_36(i,j)==0 & chl_log_215(i,j)==1)
            chl_end(i,j) = 1;
        end
    end
end

chl_first = chl_first.*tro_pac_in';
chl_both = chl_both.*tro_pac_in';
chl_end = chl_end.*tro_pac_in';

chl_first(chl_first==0) = nan;
chl_both(chl_both==0) = nan;
chl_end(chl_end==0) = nan;

chl_first(~isnan(chl_first)) = 0; 
chl_both(~isnan(chl_both)) = 1;
chl_end(~isnan(chl_end)) = 2;

C1 = [255, 128, 128]/256;
C2 = [128, 212, 255]/256;
C3 = [123, 202, 161]/256;
colortable = [C2;C3;C1];

load('D:\Geodas\coast\global_high_100more_0_360.mat')
loc=find(data_coast(:,2)==0|data_coast(:,2)==1);
num=length(loc);

R=6371;
s_fac=(2*pi*R/360*(lat(1)-lat(2)))*(2*pi*R/360*(lon(2)-lon(1))*cosd(Y));

S = (sum(sum(chl_log_215'.*tro_pac_in.*s_fac)) - sum(sum(chl_log_36'.*tro_pac_in.*s_fac)))/15
%7.16e04;
%%
figure(11)

set(gcf,'pos',[2650 250 750 400])
set(gcf,'color',[1 1 1])

hold on
plot(8:12,S_unit(1:5),'color',colortable_l(1,:),'linestyle','-','linewidth',1.15);
for i=1:19
    if sum(isnan(S_unit(6+(i-1)*12:17+(i-1)*12)))==0
    plot(1:12,S_unit(6+(i-1)*12:17+(i-1)*12),'color',colortable_l(i+1,:),'linestyle','-','linewidth',1.15);
    end
end
plot(1:7,S_unit(234:240),'color',colortable_l(21,:),'linestyle','-','linewidth',1.15);

B = datestr(time(6:17),'mmm');
D = [' 1   ';'  2  ';' 3   '];
set(gca, 'Box', 'off','layer','top','color','none', ...              % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[1 12],...
         'ylim',[1.5e7 3.5e7],...
         'Xtick',1:12,'xticklabel',B,...
         'Ytick',1e7:1e7:3e7,'yticklabel',D);
     
ylabel({'CRT area';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',12,'fontweight','normal')
xlabel('Date','FontName','Arial',...
    'fontsize',12,'fontweight','normal')
set(gca,'pos',[0.12 0.13 0.84 0.8])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*4e+07])

l_f = plot(1:12,S_first,'color',[2 38 62]/256,'linestyle','--','linewidth',1.15);
l_e = plot(1:12,S_end,'color',[2 38 62]/256,'linestyle','-','linewidth',1.15);
LG1 = legend([l_f,l_e],'average over 2003-2007','average over 2018-2022');
set(LG1,'FontName','Arial','fontsize',13,'fontweight','normal','box','off')

h1=colorbar;
colormap(colortable_l)
caxis([0 1]) 

E=datestr(time(1:12*20:end),'yyyy');

set(h1,'Ticks',[1/42:20/21:1-1/42],'TickLabels',E,'fontsize',13)
set(h1,'location','south','AxisLocation','out')

set(h1,'position',[0.16 0.22 0.5 0.03])
set(h1,'Box','on','linewidth',1.1,'TickDir', 'in', 'TickLength', [.008])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')

ax3 = axes('Position',[0.12 0.13 0.84 0.8],'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

%%

figure(12)
set(gcf,'pos',[2650 250 750 400])
set(gcf,'color',[1 1 1])

hold on

[c1,h1]=contourf(X,Y,chl_first',1,'linestyle','none');
[c2,h2]=contourf(X,Y,chl_both',1,'linestyle','none');
[c3,h3]=contourf(X,Y,chl_end',1,'linestyle','none');

h1=colorbar;
colormap(colortable)
caxis([0 2]) 

hold on
for i=1:num
    if i==num
        plot(data_coast(loc(i)+1:end,1),data_coast(loc(i)+1:end,2),'k');
    else
        plot(data_coast(loc(i)+1:loc(i+1)-1,1),data_coast(loc(i)+1:loc(i+1)-1,2),'k')
    end
    hold on
end
%对陆地进行填充
for i=1:num
    if i==num
            fix=[data_coast(loc(i)+1:end,1);data_coast(loc(i)+1,1)];
            fiy=[data_coast(loc(i)+1:end,2);data_coast(loc(i)+1,2)];
            fill(fix,fiy,[.7 .7 .7])
        else
            fix=[data_coast(loc(i)+1:loc(i+1)-1,1);data_coast(loc(i)+1,1)];
            fiy=[data_coast(loc(i)+1:loc(i+1)-1,2);data_coast(loc(i)+1,2)];
            fill(fix,fiy,[.7 .7 .7])
    end
end
axis equal
%NPSG范围
xlim([130 290])
ylim([-30 30])

B=['     ';'150°E';'180° ';'150°W';'120°W';' 90°W';'     '];
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',11,'fontweight','bold');
% D=['40°S';'    ';'20°S';'    ';' 0° ';'    ';'20°N';'    ';'40°N'];
D=['    ';'20°S';' 0° ';'20°N';'    '];
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',11,'fontweight','bold');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
grid off

set(gca,'pos',[0.08 0.24 0.88 0.6537/0.775*0.88])
set(h1,'location','southoutside');
set(h1,'Ticks',[0:2/3:2],'TickLabels',[])
% set(get(h1,'title'),'string','2003-2007 Only  overlapping of 2003-2007 and 2018-2022  2018-2022 Only',...
%     'fontsize',13,'Rotation',0,'position',[247 -30]);
set(h1,'position',[0.08 0.1439 0.88 0.04])
set(h1,'Box','on','linewidth',1.1,'TickDir', 'out', 'TickLength', [.008 .008])
% text(130+160*9/242,20,'a','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% t1=title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[133 1.03*60-30])
ax2 = axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax2,'pos',[0.08 0.29 0.88 0.6306]);

ax3 = axes('Position',[0.08 0.1439 0.88 0.04],'XAxisLocation','bottom','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'pos',[0.08 0.03 0.88 0.1]);
set(ax3,'color','none')

text(0.07,0.5,'2003-2007 Only','FontName','Arial','fontsize',13,'fontweight','normal')
text(0.35,0.5,{'          overlapping of    ';'2003-2007 and 2018-2022'},'FontName','Arial','fontsize',13,'fontweight','normal')
text(0.74,0.5,'2018-2022 Only','FontName','Arial','fontsize',13,'fontweight','normal')

set(ax3,'XTick', [],'YTick', []);   % 去掉xy轴刻度
set(ax3,'box','off');
hAxes = gca;
hAxes.YRuler.Axle.LineStyle = 'none';
hAxes.XRuler.Axle.LineStyle = 'none';

%%

figure(21)

set(gcf,'pos',[2650 250 750 680])
set(gcf,'color',[1 1 1])

s1 = subplot(2,1,1)

hold on
plot(8:12,S_unit(1:5),'color',colortable_l(1,:),'linestyle','-','linewidth',1.15);
for i=1:19
    if sum(isnan(S_unit(6+(i-1)*12:17+(i-1)*12)))==0
    plot(1:12,S_unit(6+(i-1)*12:17+(i-1)*12),'color',colortable_l(i+1,:),'linestyle','-','linewidth',1.15);
    end
end
plot(1:7,S_unit(234:240),'color',colortable_l(21,:),'linestyle','-','linewidth',1.15);

B = datestr(time(6:17),'mmm');
D = [' 1   ';'  2  ';' 3   '];
set(gca, 'Box', 'off','layer','top','color','none', ...              % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[1 12],...
         'ylim',[1.5e7 3.5e7],...
         'Xtick',1:12,'xticklabel',B,...
         'Ytick',1e7:1e7:3e7,'yticklabel',1:1:3);
     
ylabel({'CRT area';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',12,'fontweight','normal')
xlabel('Date','FontName','Arial',...
    'fontsize',12,'fontweight','normal')
set(s1,'pos',[0.11 0.6 0.85 0.34])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[1.225 3.57e+07])

% l_f = plot(1:12,S_first,'color',[2 38 62]/256,'linestyle','--','linewidth',1.15);
% l_e = plot(1:12,S_end,'color',[2 38 62]/256,'linestyle','-','linewidth',1.15);
% LG1 = legend([l_f,l_e],'average over 2003-2007','average over 2018-2022');
% set(LG1,'FontName','Arial','fontsize',13,'fontweight','normal','box','off')

h1=colorbar;
colormap(colortable_l)
caxis([0 1]) 

E=datestr(time(1:12*20:end),'yyyy');

set(h1,'Ticks',[1/42:20/21:1-1/42],'TickLabels',E,'fontsize',13)
set(h1,'location','south','AxisLocation','out')

set(h1,'position',[0.16 0.65 0.5 0.0176])
set(h1,'Box','on','linewidth',1.1,'TickDir', 'in', 'TickLength', [.008])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')

ax1 = axes('Position',[0.11 0.6 0.85 0.34],'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度


s2 = subplot(2,1,2)
hold on

[c1,h1]=contourf(X,Y,chl_first',1,'linestyle','none');
[c2,h2]=contourf(X,Y,chl_both',1,'linestyle','none');
[c3,h3]=contourf(X,Y,chl_end',1,'linestyle','none');

h2=colorbar;
colormap(colortable)
caxis([0 2]) 

hold on
for i=1:num
    if i==num
        plot(data_coast(loc(i)+1:end,1),data_coast(loc(i)+1:end,2),'k');
    else
        plot(data_coast(loc(i)+1:loc(i+1)-1,1),data_coast(loc(i)+1:loc(i+1)-1,2),'k')
    end
    hold on
end
%对陆地进行填充
for i=1:num
    if i==num
            fix=[data_coast(loc(i)+1:end,1);data_coast(loc(i)+1,1)];
            fiy=[data_coast(loc(i)+1:end,2);data_coast(loc(i)+1,2)];
            fill(fix,fiy,[.7 .7 .7])
        else
            fix=[data_coast(loc(i)+1:loc(i+1)-1,1);data_coast(loc(i)+1,1)];
            fiy=[data_coast(loc(i)+1:loc(i+1)-1,2);data_coast(loc(i)+1,2)];
            fill(fix,fiy,[.7 .7 .7])
    end
end
axis equal
%NPSG范围
xlim([130 290])
ylim([-30 30])

B_lon=['     ';'150°E';'180° ';'150°W';'120°W';' 90°W';'     '];
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B_lon,'fontsize',11,'fontweight','bold');
% D=['40°S';'    ';'20°S';'    ';' 0° ';'    ';'20°N';'    ';'40°N'];
D=['    ';'20°S';' 0° ';'20°N';'    '];
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',11,'fontweight','bold');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
grid off

set(s2,'pos',[0.11 0.11 0.85 0.3412/0.6694*0.85])
set(h2,'location','southoutside');
set(h2,'Ticks',[0:2/3:2],'TickLabels',[])
% set(get(h1,'title'),'string','2003-2007 Only  overlapping of 2003-2007 and 2018-2022  2018-2022 Only',...
%     'fontsize',13,'Rotation',0,'position',[247 -30]);
set(h2,'position',[0.11 0.075 0.85 0.0235])
set(h2,'Box','on','linewidth',1.1,'TickDir', 'out', 'TickLength', [.008 .008])
% text(130+160*9/242,20,'a','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[133 1.03*60-30])
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax2,'pos',[0.11 0.145 0.85 0.35725]);

ax3 = axes('Position',[0.11 0.03 0.85 0.02],'XAxisLocation','bottom','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'pos',[0.11 0.03 0.85 0.02]);
set(ax3,'color','none')

text(0.07,0.5,'2003-2007 Only','FontName','Arial','fontsize',13,'fontweight','normal')
text(0.34,0.5,{'          overlapping of    ';'2003-2007 and 2018-2022'},'FontName','Arial','fontsize',13,'fontweight','normal')
text(0.74,0.5,'2018-2022 Only','FontName','Arial','fontsize',13,'fontweight','normal')

set(ax3,'XTick', [],'YTick', []);   % 去掉xy轴刻度
set(ax3,'box','off');
hAxes = gca;
hAxes.YRuler.Axle.LineStyle = 'none';
hAxes.XRuler.Axle.LineStyle = 'none';

set(s1,'colormap',colortable_l);

ax4 = axes('Position',[0.6 0.7 0.3 0.2],'XAxisLocation','bottom','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax4,'pos',[0.6 0.8 0.33 0.13])

hold on
l_f = plot(1:12,S_first,'color',[2 38 62]/256,'linestyle','--','linewidth',1.15);
l_e = plot(1:12,S_end,'color',[2 38 62]/256,'linestyle','-','linewidth',1.15);

BB = datestr(time(6:17),'mmm');
BB_num = [2,3,5,6,8,9,11,12];
for i = 1:8;
    BB(BB_num(i),:) = '   ';
end
D = [' 1   ';'  2  ';' 3   '];
set(gca, 'Box', 'on','layer','top','color','none', ...              % 边框
         'linewidth', 1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[1 12],...
         'ylim',[1.8e7 3e7],...
         'Xtick',1:12,'xticklabel',BB,...
         'Ytick',2e7:1e7:3e7,'yticklabel',2:1:3);
     
ylabel({'CRT area';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal')
% xlabel('Date','FontName','Arial',...
%     'fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',10,'fontweight','normal')

LG1 = legend([l_f,l_e],'average over 2003-2007','average over 2018-2022');
set(LG1,'FontName','Arial','fontsize',10,'fontweight','normal','box','off')
set(LG1,'pos',[0.6687 0.88 0.2627 0.0537])


print(figure(21),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_4'],'-dpng','-r1200')

