
clear;clc;

filename = ['E:\data\noaa_silicate\seasonal\woa18_all_i13_01.nc'];
nitrate = ncread(filename,'i_an');
n_temp = nitrate(:,:,1);
n_temp = cat(1,n_temp(181:end,:),n_temp(1:180,:));
n_sur(:,:,1) = n_temp;

filename = ['E:\data\noaa_silicate\seasonal\woa18_all_i14_01.nc'];
nitrate = ncread(filename,'i_an');
n_temp = nitrate(:,:,1);
n_temp = cat(1,n_temp(181:end,:),n_temp(1:180,:));
n_sur(:,:,2) = n_temp;

filename = ['E:\data\noaa_silicate\seasonal\woa18_all_i15_01.nc'];
nitrate = ncread(filename,'i_an');
n_temp = nitrate(:,:,1);
n_temp = cat(1,n_temp(181:end,:),n_temp(1:180,:));
n_sur(:,:,3) = n_temp;

filename = ['E:\data\noaa_silicate\seasonal\woa18_all_i16_01.nc'];
nitrate = ncread(filename,'i_an');
n_temp = nitrate(:,:,1);
n_temp = cat(1,n_temp(181:end,:),n_temp(1:180,:));
n_sur(:,:,4) = n_temp;

lon_n = ncread(filename,'lon');
lat_n = ncread(filename,'lat');
lon_n = [lon_n(181:end);lon_n(1:180)+360];

clear nitrate n_temp filename


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

[~,month,~] = datevec(time_chl);
loc_m = find(month==12|month==1|month==2);
chl_cm(:,:,1) = nanmean(chl(:,:,loc_m),3);
loc_m = find(month==3|month==4|month==5);
chl_cm(:,:,2) = nanmean(chl(:,:,loc_m),3);
loc_m = find(month==6|month==7|month==8);
chl_cm(:,:,3) = nanmean(chl(:,:,loc_m),3);
loc_m = find(month==9|month==10|month==11);
chl_cm(:,:,4) = nanmean(chl(:,:,loc_m),3);

chl_4 = zeros(200,100,4);
for k=1:4
for i = 1:200
    for j = 1:100
        chl_4(i,j,k) = nanmean(nanmean(chl_cm(12*(i-1)+1:12*(i-1)+12,12*(j-1)+1:12*(j-1)+12,k)));
    end
end
end
lon_c1 = 100.5:1:299.5;
lat_c1 = 49.5:-1:-49.5;

clear chl_temp chl



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

u_c = u_c(:,:,116:355);
v_c = v_c(:,:,116:355);

loc_m = find(month==12|month==1|month==2);
u_cm(:,:,1) = nanmean(u_c(:,:,loc_m),3);v_cm(:,:,1) = nanmean(v_c(:,:,loc_m),3);
loc_m = find(month==3|month==4|month==5);
u_cm(:,:,2) = nanmean(u_c(:,:,loc_m),3);v_cm(:,:,2) = nanmean(v_c(:,:,loc_m),3);
loc_m = find(month==6|month==7|month==8);
u_cm(:,:,3) = nanmean(u_c(:,:,loc_m),3);v_cm(:,:,3) = nanmean(v_c(:,:,loc_m),3);
loc_m = find(month==9|month==10|month==11);
u_cm(:,:,4) = nanmean(u_c(:,:,loc_m),3);v_cm(:,:,4) = nanmean(v_c(:,:,loc_m),3);


u_4 = zeros(180,360,4);
v_4 = zeros(180,360,4);
for k=1:4;
for i = 1:180
    for j = 1:360
        u_4(i,j,k) = nanmean(nanmean(u_cm(4*(i-1)+1:4*(i-1)+3,4*(j-1)+2:4*(j-1)+4,k)));
        v_4(i,j,k) = nanmean(nanmean(v_cm(4*(i-1)+1:4*(i-1)+3,4*(j-1)+2:4*(j-1)+4,k)));
    end
end
end
lon_uv1 = 0.5:1:359.5;
lat_uv1 = -89.5:1:89.5;

clear u_c v_c


load('D:\Geodas\coast\global_high_100more_0_360.mat')
loc=find(data_coast(:,2)==0|data_coast(:,2)==1);
num=length(loc);
%%
[X_c,Y_c] = meshgrid(lon_c1,lat_c1);
colortable_c=textread('D:\Matlab\anzhuang\bin\colorbar_NCL\hotres.txt');
% colortable_c=flipud(colortable_c);
chl_temp = chl_4(:,:,4);
chl_temp(find(abs(chl_temp)>1)) = 1;

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

M1 = contour(X_c,Y_c,chl_4(:,:,1)',[0.1 0.1],'color','k','linewidth',1.2);
M2 = contour(X_c,Y_c,chl_4(:,:,2)',[0.1 0.1],'color','r','linewidth',1.2);
M3 = contour(X_c,Y_c,chl_4(:,:,3)',[0.1 0.1],'color','b','linewidth',1.2);
M4 = contour(X_c,Y_c,chl_4(:,:,4)',[0.1 0.1],'color','g','linewidth',1.2);

line_1 = M1(:,[184:340,765:895]);
line_2 = M2(:,[296:593]);
line_3 = M3(:,[479:800]);
line_4 = M4(:,[219:606]);

% plot(line_c(1,:),line_c(2,:),'color','k','linewidth',1.2);
% plot(M4(1,219:606),M4(2,219:606),'color','r','linewidth',1.2);
% plot(M1(1,184:340),M1(2,184:340),'color','r','linewidth',1.2);

%%
[X_n,Y_n] = meshgrid(double(lon_n),double(lat_n));
colortable_n=textread('D:\Matlab\anzhuang\bin\colorbar_NCL\MPL_YlOrRd.txt');
% colortable_n=flipud(colortable_n);
n_temp = n_sur;
n_temp(find(n_temp>12))=12;
[X_uv,Y_uv] = meshgrid(double(lon_uv1),double(lat_uv1));
%%
figure(21)

set(gcf,'pos',[150 160 1000 990])
set(gcf,'color',[1 1 1])

s1 = subplot(4,2,1)
n_temp = n_sur(:,:,1);
n_temp(find(n_temp>12))=12;
[c,h]=contourf(X_n,Y_n,n_temp',100,'linestyle','none');
% h2=colorbar;
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
plot(line_1(1,:),line_1(2,:),'color','k','linewidth',1.2);
B=['     ';'150°E';'180° ';'150°W';'120°W';' 90°W';'     '];
D=['    ';'20°S';' 0° ';'20°N';'    '];
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
caxis([0 8]) 
set(s1,'pos',[0.05 0.78 0.44 0.1891])
t1=title('(a) Winter','FontName','Arial','fontsize',13,'fontweight','bold','position',[145.7 1.03*60-30])
set(t1,'position',[145.7 1.03*60-30])
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax1,'pos',[0.05 0.79 0.44 0.168]);
set(s1,'FontName','Arial','fontsize',12,'fontweight','normal')

s3 = subplot(4,2,3)
n_temp = n_sur(:,:,2);
n_temp(find(n_temp>12))=12;
[c,h]=contourf(X_n,Y_n,n_temp',100,'linestyle','none');
% h2=colorbar;
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
plot(line_2(1,:),line_2(2,:),'color','k','linewidth',1.2);
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
caxis([0 8]) 
set(s3,'pos',[0.05 0.55 0.44 0.1891])
t3=title('(b) Spring','FontName','Arial','fontsize',13,'fontweight','bold','position',[135.7 1.03*60-30])
set(t3,'position',[145.7 1.03*60-30])
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax3,'pos',[0.05 0.56 0.44 0.168]);
set(s3,'FontName','Arial','fontsize',12,'fontweight','normal')

s5 = subplot(4,2,5)
n_temp = n_sur(:,:,3);
n_temp(find(n_temp>12))=12;
[c,h]=contourf(X_n,Y_n,n_temp',100,'linestyle','none');
% h2=colorbar;
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
plot(line_3(1,:),line_3(2,:),'color','k','linewidth',1.2);
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
caxis([0 8]) 
set(s5,'pos',[0.05 0.32 0.44 0.1891])
t5=title('(c) Summer','FontName','Arial','fontsize',13,'fontweight','bold','position',[135.7 1.03*60-30])
set(t5,'position',[147.7 1.03*60-30])
ax5 = axes('Position',get(s5,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax5,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax5,'pos',[0.05 0.33 0.44 0.168]);
set(s5,'FontName','Arial','fontsize',12,'fontweight','normal')

s7 = subplot(4,2,7)
n_temp = n_sur(:,:,4);
n_temp(find(n_temp>12))=12;
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
plot(line_4(1,:),line_4(2,:),'color','k','linewidth',1.2);
set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
caxis([0 8]) 
t7=title('(d) Autumn','FontName','Arial','fontsize',13,'fontweight','bold','position',[135.7 1.03*60-30])
set(t7,'position',[147.7 1.03*60-30])
set(h1,'location','southoutside');
set(h1,'Ticks',[0:2:8],'TickLabels',[0:2:8])
set(get(h1,'title'),'Rotation',0,'string','Silicate [μmol/kg]','fontsize',13,'position',[164 -30]);
set(h1,'position',[0.05 0.055 0.44 0.0125])
set(h1,'Box','on','linewidth',1.1)
set(s7,'pos',[0.05 0.09 0.44 0.1891])
ax7 = axes('Position',get(s7,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax7,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax7,'pos',[0.05 0.1 0.44 0.168]);
set(s7,'FontName','Arial','fontsize',12,'fontweight','normal')

s2 = subplot(4,2,2)
u_1 = u_4(:,:,1);v_1 = v_4(:,:,1);
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
plot(line_1(1,:),line_1(2,:),'color','k','linewidth',1.2);

set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
t2=title('(e) Winter','FontName','Arial','fontsize',13,'fontweight','bold','position',[145.7 1.03*60-30])
set(t2,'position',[145.7 1.03*60-30])
set(s2,'pos',[0.545 0.78 0.44 0.1891])%0.05 0.78 0.44 0.1891
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax2,'pos',[0.545 0.79 0.44 0.168]);
set(s2,'FontName','Arial','fontsize',12,'fontweight','normal')

s4 = subplot(4,2,4)
u_1 = u_4(:,:,2);v_1 = v_4(:,:,2);
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
plot(line_2(1,:),line_2(2,:),'color','k','linewidth',1.2);

set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
t4=title('(f) Spring','FontName','Arial','fontsize',13,'fontweight','bold','position',[145 1.03*60-30])
set(t4,'position',[145 1.03*60-30])
set(s4,'pos',[0.545 0.55 0.44 0.1891])%0.05 0.55 0.44 0.1891
ax4 = axes('Position',get(s4,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax4,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax4,'pos',[0.545 0.56 0.44 0.168]);
set(s4,'FontName','Arial','fontsize',12,'fontweight','normal')

s6 = subplot(4,2,6)
u_1 = u_4(:,:,3);v_1 = v_4(:,:,3);
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
plot(line_3(1,:),line_3(2,:),'color','k','linewidth',1.2);

set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
t6=title('(g) Summer','FontName','Arial','fontsize',13,'fontweight','bold','position',[148.1 1.03*60-30])
set(t6,'position',[148.1 1.03*60-30])
set(s6,'pos',[0.545 0.32 0.44 0.1891])
ax6 = axes('Position',get(s6,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax6,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax6,'pos',[0.545 0.33 0.44 0.168]);
set(s6,'FontName','Arial','fontsize',12,'fontweight','normal')

s8 = subplot(4,2,8)
u_1 = u_4(:,:,4);v_1 = v_4(:,:,4);
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
plot(line_4(1,:),line_4(2,:),'color','k','linewidth',1.2);

set(gca,'Xtick',[130,150:30:270,290],'xticklabel',B,'fontsize',12,'fontweight','normal');
set(gca,'Ytick',[-30,-20,0,20,30],'Yticklabel',D,'fontsize',12,'fontweight','normal');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
t8=title('(h) Autumn','FontName','Arial','fontsize',13,'fontweight','bold','position',[148.1 1.03*60-30])
set(t8,'position',[148.1 1.03*60-30])
set(s8,'pos',[0.545 0.09 0.44 0.1891])
ax8 = axes('Position',get(s8,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax8,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax8,'pos',[0.545 0.1 0.44 0.168])
set(s8,'FontName','Arial','fontsize',12,'fontweight','normal')

axe = axes('Position',[0.545 0.01 0.44 0.1891/2],'XAxisLocation','bottom','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(axe,'pos',[0.545 0.01 0.44 0.1891/2]);
set(axe,'color','none')
Qc=quiver(248.5,16,0.825,0,11,'filled');
xlim([130 290])
ylim([0 30])
Qc.LineWidth = 1.2;
Qc.MaxHeadSize = 1;
Qc.AutoScale = 'on';
Qc.AutoScaleFactor = 11;
text(225,16,'current','FontName','Arial','fontsize',14,'fontweight','normal')
text(260,16,'0.3 m/s','FontName','Arial','fontsize',14,'fontweight','normal')
set(axe,'XTick', [],'YTick', []);   % 去掉xy轴刻度
set(axe,'box','off');
set(axe,'color','none')
hAxes = gca;
hAxes.YRuler.Axle.LineStyle = 'none';
hAxes.XRuler.Axle.LineStyle = 'none';


print(figure(21),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_2'],'-dpng','-r1200')

