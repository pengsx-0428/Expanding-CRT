
clear;clc;

filename = ['E:\data\noaa_silicate\annual\woa18_all_i00_01.nc'];
lon_n = ncread(filename,'lon');
lat_n = ncread(filename,'lat');
lon_n = [lon_n(181:end);lon_n(1:180)+360];
nitrate = ncread(filename,'i_an');
n_sur = nitrate(:,:,1);
n_sur = cat(1,n_sur(181:end,:),n_sur(1:180,:));

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
lon_c=ncread([filepath_chl,playname{1}],'lon');
lon_c=cat(1,lon_c(2161:end),lon_c(1:2160)+360);
lon_c=lon_c(1201:3600);
%lat:50S-50N:481:1680
lat_c=ncread([filepath_chl,playname{1}],'lat');
lat_c=lat_c(481:1680);

chl_cm = nanmean(chl,3);

chl_1 = zeros(200,100);

for i = 1:200
    for j = 1:100
        chl_1(i,j) = nanmean(nanmean(chl_cm(12*(i-1)+1:12*(i-1)+12,12*(j-1)+1:12*(j-1)+12)));
    end
end
lon_c1 = 100.5:1:299.5;
lat_c1 = 49.5:-1:-49.5;

file_uv = ['E:\data\OSCAR\data_mm_1_4_global\'];
dir_uv = dir(fullfile(file_uv,'*.mat'));
playname_uv = {dir_uv.name};
for i=1:length(playname_uv)
    year_uv(i,:) = playname_uv{i}(19:22);
end

for i=1:length(playname_uv)
    load([file_uv,playname_uv{i}]);
    
        u_temp = u;
        v_temp = v;
    if i==1
        u_c=u_temp;
        v_c=v_temp;
    else
        u_c = cat(3,u_c,u_temp);
        v_c = cat(3,v_c,v_temp);
    end
end
clear u_temp v_temp ug vg u v
lon_uv = ncread('E:\data\OSCAR\oscar_currents_final_19930101.nc','lon');
lat_uv = ncread('E:\data\OSCAR\oscar_currents_final_19930101.nc','lat');
for i=1:length(playname_uv)
    for j = 1:12
        time_uv((i-1)*12+j) = datenum([year_uv(i,:),sprintf('%02d',j),'01'],'yyyymmdd');
    end
end

u_cm = nanmean(u_c(:,:,116:355),3);
v_cm = nanmean(v_c(:,:,116:355),3);

u_1 = zeros(180,360);
v_1 = zeros(180,360);

for i = 1:180
    for j = 1:360
        u_1(i,j) = nanmean(nanmean(u_cm(4*(i-1)+1:4*(i-1)+3,4*(j-1)+2:4*(j-1)+4)));
        v_1(i,j) = nanmean(nanmean(v_cm(4*(i-1)+1:4*(i-1)+3,4*(j-1)+2:4*(j-1)+4)));
    end
end
lon_uv1 = 0.5:1:359.5;
lat_uv1 = -89.5:1:89.5;

load('D:\Geodas\coast\global_high_100more_0_360.mat')
loc=find(data_coast(:,2)==0 | data_coast(:,2)==1);
num=length(loc);

[X_c,Y_c] = meshgrid(lon_c1,lat_c1);
colortable_c=textread('D:\Matlab\anzhuang\bin\colorbar_NCL\hotres.txt');
chl_temp = chl_1;
chl_temp(find(abs(chl_temp)>1)) = 1;

[X_n,Y_n] = meshgrid(double(lon_n),double(lat_n));
colortable_n=textread('D:\Matlab\anzhuang\bin\colorbar_NCL\MPL_YlOrRd.txt');
n_temp = n_sur;
n_temp(find(n_temp>12))=12;

[X_uv,Y_uv] = meshgrid(double(lon_uv1),double(lat_uv1));


%%

figure(11)
set(gcf,'pos',[2650 400 800 400])
set(gcf,'color',[1 1 1])

[c,h]=contourf(X_c,Y_c,chl_temp',500,'linestyle','none');
h1=colorbar;
colormap(colortable_c)
hold on
for i=1:num
    if i==num
        plot(data_coast(loc(i)+1:end,1),data_coast(loc(i)+1:end,2),'k');
    else
        plot(data_coast(loc(i)+1:loc(i+1)-1,1),data_coast(loc(i)+1:loc(i+1)-1,2),'k')
    end
    hold on
end
axis equal
%NPSG范围
xlim([130 290])
ylim([-30 30])

M = contour(X_c,Y_c,chl_temp',[0.1 0.1],'color','k','linewidth',1.2);
%
line_c = M(:,232:518);
plot(line_c(1,:),line_c(2,:),'color','r','linewidth',1.2);

%%
figure(12)
set(gcf,'pos',[2650 400 800 400])
set(gcf,'color',[1 1 1])

[c,h]=contourf(X_n,Y_n,n_temp',100,'linestyle','none');
h1=colorbar;
colormap(colortable_n)
hold on
for i=1:num
    if i==num
        plot(data_coast(loc(i)+1:end,1),data_coast(loc(i)+1:end,2),'k');
    else
        plot(data_coast(loc(i)+1:loc(i+1)-1,1),data_coast(loc(i)+1:loc(i+1)-1,2),'k')
    end
    hold on
end
axis equal
%NPSG范围
xlim([130 290])
ylim([-30 30])

caxis([0 7]) 

plot(line_c(1,:),line_c(2,:),'color','k','linewidth',1.2);

%%

figure(21)

set(gcf,'pos',[2650 200 700 850])
set(gcf,'color',[1 1 1])

s1 = subplot(3,1,1)
[c,h]=contourf(X_c,Y_c,chl_temp',500,'linestyle','none');
h1=colorbar;
colormap(colortable_c)
hold on
for i=1:num
    if i==num
        plot(data_coast(loc(i)+1:end,1),data_coast(loc(i)+1:end,2),'k');
    else
        plot(data_coast(loc(i)+1:loc(i+1)-1,1),data_coast(loc(i)+1:loc(i+1)-1,2),'k')
    end
    hold on
end
axis equal
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
%范围
xlim([130 290])
ylim([-30 30])
plot(line_c(1,:),line_c(2,:),'color','k','linewidth',1.2);
% B=['     ';'120°E';'     ';'140°E';'     ';'160°E';'     ';'180°E';...
%     '     ';'160°W';'     ';'140°W';'     ';'120°W';'     ';'100°W';'     ';' 80°W';'     ';];
B=['     ';'150°E';'180° ';'150°W';'120°W';' 90°W';'     '];
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',11,'fontweight','bold');
% D=['40°S';'    ';'20°S';'    ';' 0° ';'    ';'20°N';'    ';'40°N'];
D=['    ';'20°S';' 0° ';'20°N';'    '];
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',11,'fontweight','bold');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
grid off

caxis([0 1]) 
% title('Climatological Jan map of current','fontsize',14,'fontweight','bold')% set(h1,'location','southoutside')
set(h1,'Ticks',[0:0.2:1],'TickLabels',[0:0.2:1],'linewidth',1.1)
set(get(h1,'title'),'Rotation',90,'string','Chl [mg/m^3]','fontsize',13,'position',[42 80]);
set(h1,'position',[0.91 0.7048 0.0125 0.2563])
set(gca,'pos',[0.075 0.7 0.83 0.2157/0.6735*0.83])
% text(130+160*9/242,20,'a','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t1=title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[133 1.03*60-30])
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax1,'pos',[0.075 0.706 0.83 0.255]);

s2 = subplot(3,1,2)
Q=quiver(X_uv(1:4:end,1:4:end),Y_uv(1:4:end,1:4:end),u_1(1:4:end,1:4:end),v_1(1:4:end,1:4:end),11);

Q.MaxHeadSize = 0.01;
Q.AutoScale = 'on';
Q.AutoScaleFactor = 3;

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
plot(line_c(1,:),line_c(2,:),'color','k','linewidth',1.2);

set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',11,'fontweight','bold');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',11,'fontweight','bold');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
grid off
% title('Climatological Jan map of current','fontsize',14,'fontweight','bold')% set(h1,'location','southoutside')
set(gca,'pos',[0.075 0.375 0.83 0.2157/0.6735*0.83])
% text(130+160*9/242,20,'c','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% ax1_pos = ;
% ax1_pos = ax1_pos + [0.05 -0.05 -0.05 -0.05];
t3=title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[133 1.03*60-30])
set(t3,'pos',[133 1.03*60-30])
% text(230,-1,'SEC','color','r','FontName','Arial','fontsize',15,'fontweight','bold')
line([170,260,260,170,170],[-5,-5,5,5,-5],'color',[249,41,42]/256,'linewidth',1.2,'linestyle','--')
ax3 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax3,'pos',[0.075 0.053+0.325 0.83 0.258]);

s3 = subplot(3,1,3)
[c,h]=contourf(X_n,Y_n,n_temp',100,'linestyle','none');
h2=colorbar;
colormap(colortable_n)
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
%NPSG范围
axis equal
%NPSG范围
xlim([130 290])
ylim([-30 30])
plot(line_c(1,:),line_c(2,:),'color','k','linewidth',1.2);
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',11,'fontweight','bold');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',11,'fontweight','bold');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
grid off
caxis([0 8]) 
% title('Climatological Jan map of current','fontsize',14,'fontweight','bold')% set(h1,'location','southoutside')
set(h2,'Ticks',[0:2:12],'TickLabels',[0:2:12],'linewidth',1.1)
set(get(h2,'title'),'Rotation',90,'string','Silicate [μmol/kg]','fontsize',13,'position',[38 80]);
set(h2,'position',[0.91 0.0548 0.0125 0.2563])
set(gca,'pos',[0.075 0.05 0.83 0.2157/0.6735*0.83])
% text(130+160*9/242,20,'b','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t2=title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[133 1.03*60-30])
ax2 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax2,'pos',[0.075 0.38-0.325 0.83 0.256]);

ax4 = axes('Position',[0.6 0.3+0.325 0.3 0.2],'XAxisLocation','bottom','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax4,'pos',[0.6 0.312+0.325 0.35 0.05]);
set(ax4,'color','none')

Qc=quiver(30,4,0.5,0,11,'filled');
xlim([0 160/0.83*0.312])
ylim([0 60/0.258*0.05])
Qc.LineWidth = 1.2;
Qc.MaxHeadSize = 1;
Qc.AutoScale = 'on';
Qc.AutoScaleFactor = 11;
text(38,4,'0.3 m/s','FontName','Arial','fontsize',13,'fontweight','normal')
text(14,4,'current','FontName','Arial','fontsize',13,'fontweight','normal')
set(ax4,'XTick', [],'YTick', []);   % 去掉xy轴刻度
set(ax4,'box','off');
hAxes = gca;
hAxes.YRuler.Axle.LineStyle = 'none';
hAxes.XRuler.Axle.LineStyle = 'none';

set(s1,'colormap',colortable_c);

print(figure(21),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_1'],'-dpng','-r1200')




