
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

% load('E:\matlab_code\Tropical_Pacific\tropical_high_chl_region2.mat');
% highc_re_n=[line_west;flipud(line_east);line_west(1,:)];
highc_re_n=[120,120,260,260,120;5,40,40,5,5]';
highc_re_s=[160,160,280,280,160;-40,-5,-5,-40,-40]';
% highc_re_s=[160,160,240,240,160;-40,-5,-5,-40,-40]';

load('E:\data\chl_9km_monthly\occci_in.mat');
in_temp=cat(2,in_temp(:,2161:end),in_temp(:,1:2160));
sea_in=~in_temp(481:1680,1201:3600);
x=x(481:1680,1201:3600);y=y(481:1680,1201:3600);

%%
tro_pac_in_n=inpolygon(X,Y,highc_re_n(:,1),highc_re_n(:,2));
tro_pac_in_s=inpolygon(X,Y,highc_re_s(:,1),highc_re_s(:,2));
load('E:\matlab_code\Tropical_Pacific\tropical_high_chl_region2.mat');
highc_re_in=[line_west;flipud(line_east);line_west(1,:)];
tro_pac_in_eq=inpolygon(X,Y,highc_re_in(:,1),highc_re_in(:,2));

for i=1:240
    in_chl07=zeros(2400,1200);
    chl_temp=chl(:,:,i);
    for j=1:2400
        for k=1:1200
            if chl_temp(j,k)<0.07;
                in_chl07(j,k)=1;
            end
        end
    end
    inploy_temp = sea_in.*tro_pac_in_n.*in_chl07';
%     inploy_temp(inploy_temp==0) =nan;
    num_temp =2400*1200 - length(find(inploy_temp==0));
    lon_npsg_m(i)=sum(sum(X.*inploy_temp))/num_temp;
    lat_npsg_m(i)=sum(sum(Y.*inploy_temp))/num_temp;
end

for i=1:240
    in_chl07=zeros(2400,1200);
    chl_temp=chl(:,:,i);
    for j=1:2400
        for k=1:1200
            if chl_temp(j,k)<0.07;
                in_chl07(j,k)=1;
            end
        end
    end
    inploy_temp = sea_in.*tro_pac_in_s.*in_chl07';
%     inploy_temp(inploy_temp==0) =nan;
    num_temp =2400*1200 - length(find(inploy_temp==0));
    lon_spsg_m(i)=sum(sum(X.*inploy_temp))/num_temp;
    lat_spsg_m(i)=sum(sum(Y.*inploy_temp))/num_temp;
end

for i=1:240
    in_chl01=zeros(2400,1200);
    chl_temp=chl(:,:,i);
    for j=1:2400
        for k=1:1200
            if chl_temp(j,k)>0.1;
                in_chl01(j,k)=1;
            end
        end
    end
    inploy_temp = sea_in.*tro_pac_in_eq.*in_chl01';
    num_temp =2400*1200 - length(find(inploy_temp==0));
    lon_eq_m(i)=sum(sum(X.*inploy_temp))/num_temp;
    lat_eq_m(i)=sum(sum(Y.*inploy_temp))/num_temp;
end

[imf_2, residual_2] = emd(lon_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lon_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lon_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lon_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lon_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',6);
Tt_no(:,1) = residual_2;
Tt_no(:,2) = residual_3;
Tt_no(:,3) = residual_4;
Tt_no(:,4) = residual_5;
Tt_no(:,5) = residual_6;
for i=1:240
    lon_npsg_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    lon_npsg_m_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end

[imf_2, residual_2] = emd(lat_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lat_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lat_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lat_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lat_npsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',6);
Tt_na(:,1) = residual_2;
Tt_na(:,2) = residual_3;
Tt_na(:,3) = residual_4;
Tt_na(:,4) = residual_5;
Tt_na(:,5) = residual_6;
for i=1:240
    lat_npsg_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    lat_npsg_m_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end

Tt_n_lon_yr = (lon_npsg_m_m(end) - lon_npsg_m_m(1))/240*12
Tt_n_lon_m_2 = lon_npsg_m_m - detrend(lon_npsg_m_m);
Tt_n_lon_yr_2 = (Tt_n_lon_m_2(end) - Tt_n_lon_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_no(:,i) - detrend(Tt_no(:,i));
    Tt_no_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_no = std(Tt_no_all);

Tt_n_lat_yr = (lat_npsg_m_m(end) - lat_npsg_m_m(1))/240*12
Tt_n_lat_m_2 = lat_npsg_m_m - detrend(lat_npsg_m_m);
Tt_n_lat_yr_2 = (Tt_n_lat_m_2(end) - Tt_n_lat_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_na(:,i) - detrend(Tt_na(:,i));
    Tt_na_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_na = std(Tt_na_all);

[imf_2, residual_2] = emd(lon_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lon_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lon_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lon_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lon_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',6);
Tt_so(:,1) = residual_2;
Tt_so(:,2) = residual_3;
Tt_so(:,3) = residual_4;
Tt_so(:,4) = residual_5;
Tt_so(:,5) = residual_6;
for i=1:240
    lon_spsg_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    lon_spsg_m_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end

[imf_2, residual_2] = emd(lat_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(lat_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(lat_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(lat_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(lat_spsg_m,'SiftRelativeTolerance',0,'SiftMaxIterations',6);
Tt_sa(:,1) = residual_2;
Tt_sa(:,2) = residual_3;
Tt_sa(:,3) = residual_4;
Tt_sa(:,4) = residual_5;
Tt_sa(:,5) = residual_6;
for i=1:240
    lat_spsg_std(i,1) = std([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end
for i=1:240
    lat_spsg_m_m(i,1) = nanmean([residual_5(i),residual_6(i),residual_2(i),residual_3(i),...
        residual_4(i)]);
end

Tt_s_lon_yr = (lon_spsg_m_m(end) - lon_spsg_m_m(1))/240*12
Tt_s_lon_m_2 = lon_spsg_m_m - detrend(lon_spsg_m_m);
Tt_s_lon_yr_2 = (Tt_s_lon_m_2(end) - Tt_s_lon_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_so(:,i) - detrend(Tt_so(:,i));
    Tt_so_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_so = std(Tt_so_all);

Tt_s_lat_yr = (lat_spsg_m_m(end) - lat_spsg_m_m(1))/240*12
Tt_s_lat_m_2 = lat_spsg_m_m - detrend(lat_spsg_m_m);
Tt_s_lat_yr_2 = (Tt_s_lat_m_2(end) - Tt_s_lat_m_2(1))/240*12
for i=1:5
    Tt_dt = Tt_sa(:,i) - detrend(Tt_sa(:,i));
    Tt_sa_all(i) = (Tt_dt(end) - Tt_dt(1))/240*12;
end
Tt_std_sa = std(Tt_sa_all);

%%
chl_cm = nanmean(chl,3);
chl_1 = zeros(200,100);

for i = 1:200
    for j = 1:100
        chl_1(i,j) = nanmean(nanmean(chl_cm(12*(i-1)+1:12*(i-1)+12,12*(j-1)+1:12*(j-1)+12)));
    end
end
lon_c1 = 100.5:1:299.5;
lat_c1 = 49.5:-1:-49.5;

load('D:\Geodas\coast\global_high_100more_0_360.mat')
loc=find(data_coast(:,2)==0|data_coast(:,2)==1);
num=length(loc);
% C1 = [255,31,91]/256;
% C2 = [0,154,222]/256;
[X_c,Y_c] = meshgrid(lon_c1,lat_c1);
colortable_c=textread('D:\Matlab\anzhuang\bin\colorbar_NCL\hotres.txt');

chl_temp = chl_1;
chl_temp(find(abs(chl_temp)>1)) = 1;

M1 = contour(X_c,Y_c,chl_temp',[0.1 0.1],'color','k','linewidth',1.2);
M7 = contour(X_c,Y_c,chl_temp',[0.07 0.07],'color','r','linewidth',1.2);

%
line_eq = M1(:,248:497);
plot(line_eq(1,:),line_eq(2,:),'color','r','linewidth',1.2);
line_n = M7(:,68:380);
plot(line_n(1,:),line_n(2,:),'color','k','linewidth',1.2);
line_s = M7(:,482:816);
plot(line_s(1,:),line_s(2,:),'color','k','linewidth',1.2);
caxis([0 1]) 

%%
    in_chl07=zeros(2400,1200);
    chl_temp=chl_cm;
    for j=1:2400
        for k=1:1200
            if chl_cm(j,k)<0.07;
                in_chl07(j,k)=1;
            end
        end
    end
    inploy_s_temp = sea_in.*tro_pac_in_s.*in_chl07';
    inploy_n_temp = sea_in.*tro_pac_in_n.*in_chl07';
%     inploy_temp(inploy_temp==0) =nan;
    num_n_temp =2400*1200 - length(find(inploy_n_temp==0));
    num_s_temp =2400*1200 - length(find(inploy_s_temp==0));
    c_npsg(1)=sum(sum(X.*inploy_n_temp))/num_n_temp;
    c_npsg(2)=sum(sum(Y.*inploy_n_temp))/num_n_temp;
    
    c_spsg(1)=sum(sum(X.*inploy_s_temp))/num_s_temp;
    c_spsg(2)=sum(sum(Y.*inploy_s_temp))/num_s_temp;

C1=[0, 0, 200]/256;
C2=[0 176 0]/256;
C3=[233 0 45]/256;

%%


scatter(c_npsg(1),c_npsg(2),50,'filled','MarkerFaceColor',C2)
Qn=quiver(c_npsg(1),c_npsg(2),lon_npsg_m_m(end)-lon_npsg_m_m(1),lat_npsg_m_m(end)-lat_npsg_m_m(1),11);
Qn.LineWidth = 3;
Qn.MaxHeadSize = 0.6;
Qn.AutoScale = 'on';
Qn.AutoScaleFactor = 11;
Qn.Color = C2;
scatter(c_spsg(1),c_spsg(2),50,'filled','MarkerFaceColor',C2)
Qs=quiver(c_spsg(1),c_spsg(2),lon_spsg_m_m(end)-lon_spsg_m_m(1),lat_spsg_m_m(end)-lat_spsg_m_m(1),11);
Qs.LineWidth = 3;
Qs.MaxHeadSize = 0.6;
Qs.AutoScale = 'on';
Qs.AutoScaleFactor = 11;
Qs.Color = C2;
%%

figure(21)
set(gcf,'pos',[2650 200 800 840])
set(gcf,'color',[1 1 1])

s1=subplot(3,2,1)
hold on
fn1 = fill([time_chl';flipud(time_chl')],[lon_npsg_m_m+2*lon_npsg_std;flipud(lon_npsg_m_m-2*lon_npsg_std)],[.7 .7 .7])
set(fn1,'EdgeColor','none','FaceAlpha',0.5)
ln2 = plot(time_chl,lon_npsg_m,'color',[2 38 62]/256,'linewidth',1.3)
ln1 = plot(time_chl,lon_npsg_m_m,'color',[233 0 45]/256,'linewidth',1.2);
B=datestr(time_chl(6+12:48:end),'yyyy');
C = ['165°E';'180° ';'165°W'];
set(gca, 'Box', 'off','layer','top','color','none', ...         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[165 195],...
         'Xtick',time_chl(18:48:end),'xticklabel',B,...
         'Ytick',165:15:195,'yticklabel',C);
ylabel('Longitude','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t1 = title('(a) Longitude of NPSG central location','fontsize',13,'fontweight','bold');
set(t1,'pos',[time_chl(1)+3900 1.03*30+165])
set(s1,'pos',[0.11 0.74 0.37 0.21])
% t1 = text(time_chl(150),170,[sprintf('%5.3f',Tt_n_lon_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t2 = text(time_chl(96),167,[sprintf('%5.2f',Tt_n_lon_yr_2*240/12),'(±',sprintf('%4.2f',Tt_std_no*240/12),') °, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',12,'fontweight','normal');
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s2=subplot(3,2,2)
hold on
fn2 = fill([time_chl';flipud(time_chl')],[lat_npsg_m_m+2*lat_npsg_std;flipud(lat_npsg_m_m-2*lat_npsg_std)],[.7 .7 .7])
set(fn2,'EdgeColor','none','FaceAlpha',0.5)
ln4 = plot(time_chl,lat_npsg_m,'color',[2 38 62]/256,'linewidth',1.3)
ln3 = plot(time_chl,lat_npsg_m_m,'color',[233 0 45]/256,'linewidth',1.2);

D = ['14°N';'19°N';'24°N'];
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[14 24],...
         'Xtick',time_chl(18:48:end),'xticklabel',B,...
         'Ytick',14:5:24,'yticklabel',D);
ylabel('Latitude','FontName','Arial','fontsize',12,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t2 = title('(b) Latitude of NPSG central location','fontsize',13,'fontweight','bold');
set(t2,'pos',[time_chl(1)+3690 1.03*10+14])
set(s2,'pos',[0.6 0.74 0.37 0.21])
% t3 = text(time_chl(150),15.2,[sprintf('%5.3f',Tt_n_lat_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time_chl(96),14+10/15,[sprintf('%5.2f',Tt_n_lat_yr_2*240/12),'(±',sprintf('%4.2f',Tt_std_na*240/12),') °, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',12,'fontweight','normal');
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s3=subplot(3,2,3)
hold on
fs1 = fill([time_chl';flipud(time_chl')],[lon_spsg_m_m+2*lon_spsg_std;flipud(lon_spsg_m_m-2*lon_spsg_std)],[.7 .7 .7])
set(fs1,'EdgeColor','none','FaceAlpha',0.5)
ls2 = plot(time_chl,lon_spsg_m,'color',[2 38 62]/256,'linewidth',1.3)
ls1 = plot(time_chl,lon_spsg_m_m,'color',[233 0 45]/256,'linewidth',1.2);

E = ['150°W';'140°W';'130°W'];
set(gca, 'Box', 'off','layer','top','color','none', ...                 % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[210 230],...
         'Xtick',time_chl(18:48:end),'xticklabel',B,...
         'Ytick',210:10:230,'yticklabel',E);

ylabel('Longitude','FontName','Arial','fontsize',12,'fontweight','bold')
set(s3,'FontName','Arial','fontsize',12,'fontweight','normal')
t3 = title('(c) Longitude of SPSG central location','fontsize',13,'fontweight','bold');
set(t3,'pos',[time_chl(1)+3880 1.03*20+210])
set(s3,'pos',[0.11 0.455 0.37 0.21])
% t5 = text(time_chl(150+12),228,[sprintf('%5.3f',Tt_s_lon_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t6 = text(time_chl(96),210+20/15,[sprintf('%5.2f',Tt_s_lon_yr_2*240/12),'(±',sprintf('%4.2f',Tt_std_so*240/12),') °, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',12,'fontweight','normal');
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s4=subplot(3,2,4)
hold on
fs2 = fill([time_chl';flipud(time_chl')],[lat_spsg_m_m+2*lat_spsg_std;flipud(lat_spsg_m_m-2*lat_spsg_std)],[.7 .7 .7])
set(fs2,'EdgeColor','none','FaceAlpha',0.5)
ls4 = plot(time_chl,lat_spsg_m,'color',[2 38 62]/256,'linewidth',1.3)
ls3 = plot(time_chl,lat_spsg_m_m,'color',[233 0 45]/256,'linewidth',1.2);

F = ['28°S';'23°S';'18°S'];
set(gca, 'Box', 'off','layer','top','color','none', ...            % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time_chl(1) time_chl(end)],...
         'ylim',[-28 -18],...
         'Xtick',time_chl(18:48:end),'xticklabel',B,...
         'Ytick',-28:5:-18,'yticklabel',F);

ylabel('Latitude','FontName','Arial','fontsize',12,'fontweight','bold')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t4 = title('(d) Latitude of SPSG central location','fontsize',13,'fontweight','bold');
set(t4,'pos',[time_chl(1)+3700 1.03*10-28])
set(s4,'pos',[0.6 0.455 0.37 0.21])
% t7 = text(time_chl(150),-26,[sprintf('%5.3f',Tt_s_lat_yr),' °/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t8 = text(time_chl(96),10/15-28,[sprintf('%5.2f',Tt_s_lat_yr_2*240/12),'(±',sprintf('%4.2f',Tt_std_sa*240/12),') °, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',12,'fontweight','normal');
ax4 = axes('Position',get(s4,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax4,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
%%
chl_temp = chl_1;
chl_temp(find(chl_temp>1))=1;
s5=subplot(3,2,5)
hold on
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
%NPSG范围
xlim([110 290])
ylim([-40 40])

caxis([0 1]) 
B=['     ';'120°E';'150°E';'180° ';'150°W';'120°W';' 90°W';'     '];
set(gca,'Xtick',[110,120:30:270,290],'xticklabel',B,'fontsize',11,'fontweight','bold');
D=['40°S';'20°S';' 0° ';'20°N';'40°N'];
set(gca,'Ytick',[-40:20:40],'Yticklabel',D,'fontsize',11,'fontweight','bold');
set(gca,'layer','top', 'Box', 'off','linewidth', 1.2,'TickDir', 'out', 'TickLength', [.01 .01])
C1=[0, 83, 233]/256;
p1 = plot(line_n(1,:),line_n(2,:),'color',C1,'linewidth',1.5);
p2 = plot(line_s(1,:),line_s(2,:),'color',C2,'linewidth',1.5);
p3 = plot(line_eq(1,:),line_eq(2,:),'color',C3,'linewidth',1.5);

set(s5,'pos',[0.11 -0.15 0.79 0.7288])
t5 = title('(e)','fontsize',13,'fontweight','bold');
set(t5,'pos',[113 1.03*80-40])
set(h1,'Ticks',[0:0.2:1],'TickLabels',[0:0.2:1],'linewidth',1.1)
set(get(h1,'title'),'Rotation',90,'string','Chl [mg/m^3]','fontsize',13,'fontweight','normal','position',[44 108]);
set(h1,'position',[0.91 0.0471 0.0125 0.3344])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')


scatter(c_npsg(1),c_npsg(2),45,'filled','MarkerFaceColor',C1)
Qn=quiver(c_npsg(1),c_npsg(2),lon_npsg_m_m(end)-lon_npsg_m_m(1),lat_npsg_m_m(end)-lat_npsg_m_m(1),8);
Qn.LineWidth = 3;
Qn.MaxHeadSize = 0.6;
Qn.AutoScale = 'on';
Qn.AutoScaleFactor = 8;
Qn.Color = C1;
scatter(c_spsg(1),c_spsg(2),45,'filled','MarkerFaceColor',C2)
Qs=quiver(c_spsg(1),c_spsg(2),lon_spsg_m_m(end)-lon_spsg_m_m(1),lat_spsg_m_m(end)-lat_spsg_m_m(1),8);
Qs.LineWidth = 3;
Qs.MaxHeadSize = 0.6;
Qs.AutoScale = 'on';
Qs.AutoScaleFactor = 8;
Qs.Color = C2;
Qe=quiver(220,0,-7.2,0,8);
Qe.LineWidth = 3;
Qe.MaxHeadSize = 0.6;
Qe.AutoScale = 'on';
Qe.AutoScaleFactor = 4;
Qe.Color = C3;

ax1 = axes('Position',get(s5,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.2);   % 去掉xy轴刻度
set(ax1,'pos',[0.11 0.0468 0.79 0.335]);

print(figure(21),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_5'],'-dpng','-r1200')

