
clear;clc;

filename =['E:\data\ENSO\PDO.txt']; %grd 文件名

fid=fopen(filename,'r');
for i=1:171
    str2=fgetl(fid);
    A(i,:)=strread(str2,'%f','delimiter',' ');
end
fclose(fid);

for i=1:171
    for j=1:12
        time_pdo((i-1)*12+j)=datenum([sprintf('%4d',A(i,1)),'-',sprintf('%2d',j)],'yyyy-mm');
    end
end

pdo = reshape(A(:,2:end)',[171*12,1]);
time_pdo(pdo>50) = [];
time_pdo = time_pdo';
pdo(pdo>50) = [];

load('E:\进展\文章相关\EEP_HCT\nc_reviewer\code\S_crt.mat');
S_ext = Extend_sig_v2(S_crt(1:240),{'asymw'},600);

%%
fs = 1;
fc_low = 1/(8*12);
N = 1;
[b,a] = butter(N,fc_low/(fs/2),'low');

pdo_8 = filter(b,a,pdo);
S_8 = filter(b,a,S_crt);
S_e_8 = filter(b,a,S_ext);
S_e_8 = S_e_8(606:847);

C1 = [255, 128, 128]/256;
C2 = [128, 212, 255]/256;
pdo_p = pdo_8;pdo_p(pdo_8<0)=0;
pdo_n = pdo_8;pdo_n(pdo_8>0)=0;

%%

[imf_2, residual_2] = emd(S_ext,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(S_ext,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(S_ext,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(S_ext,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(S_ext,'SiftRelativeTolerance',0,'SiftMaxIterations',6);


It(:,1) = imf_2(:,3) + imf_2(:,4);
It(:,2) = imf_3(:,3) + imf_3(:,4);
It(:,3) = imf_4(:,3) + imf_4(:,4);
It(:,4) = imf_5(:,3) + imf_5(:,4);
It(:,5) = imf_6(:,3) + imf_6(:,4);
It_m = mean(It,2);
It_m = It_m(606:847,:);

MA(:,1) = imf_2(:,5) + imf_2(:,6) + imf_2(:,7);
MA(:,2) = imf_3(:,5) + imf_3(:,6) + imf_3(:,7);
MA(:,3) = imf_4(:,5) + imf_4(:,6) + imf_4(:,7);
MA(:,4) = imf_5(:,5) + imf_5(:,6) + imf_5(:,7);
MA(:,5) = imf_6(:,5) + imf_6(:,6) + imf_6(:,7) + imf_6(:,8);

Tt = [residual_2 , residual_3 , residual_4 , residual_5 , residual_6];
Tt = Tt(601:840,:);
Tt_m = mean(Tt,2);
for i=1:240
Tt_std(i) = std(Tt(i,:));
end

Tt_m_2 = Tt_m - detrend(Tt_m);
Tt_yr_2 = (Tt_m_2(end) - Tt_m_2(1))/240*12
for i=1:5
    Tt_m_s = Tt(:,i) - detrend(Tt(:,i));
    Tt_m_all(i) = (Tt_m_s(end) - Tt_m_s(1))/240*12;
end
Tt_s_std = std(Tt_m_all);

for i=1:1440
    MA_m(i) = nanmean(MA(i,:));
    MA_std(i) = std(MA(i,:));
end

MA_m = MA_m(601:840);
MA_std = MA_std(601:840);

[r_MA,p_MA] = corrcoef(MA_m,pdo_8(1784:2023))
% [r_LOW,p_LOW] = corrcoef(S_e_8,pdo_8(1784:2023))

%%
imf_p = imf_6;
period_imf = nan(size(imf_p,2),1);

for i = 1:size(imf_p,2)
    pan = imf_p(1:end-1,i).*imf_p(2:end,i);
    loc = find(pan<0);
    if length(loc) == 0;
    else
        zero_loc =loc + (abs(imf_p(loc,i)))./(abs(imf_p(loc,i))+abs(imf_p(loc+1,i)));
        period_imf(i) = nanmean(zero_loc(2:end) - zero_loc(1:end-1))/12*2;
    end
end

%%

figure(11)

set(gcf,'pos',[2650 200 770 500])
set(gcf,'color',[1 1 1])

s1=subplot(2,1,1)
hold on
yyaxis right
a1 = area(time_pdo,pdo_p,'LineWidth',1,'FaceColor',C1,'EdgeColor','none','FaceAlpha',.6,'EdgeAlpha',1,'ShowBaseLine','off');
hold on
a2 = area(time_pdo,pdo_n,'LineWidth',1,'FaceColor',C2,'EdgeColor','none','FaceAlpha',.6,'EdgeAlpha',1,'ShowBaseLine','off');
xlim([time(1) time(end-2)])
ylim([-3 3])

set(gca,'Xtick',time(18:48:end),'xticklabel',[],'fontsize',12,'fontweight','normal');
set(gca,'Ytick',-3:3:3,'yticklabel',-3:3:3,'fontsize',12,'fontweight','normal');
ylabel('PDO index','Rotation',90,'FontName','Arial','fontsize',11,...
    'fontweight','normal','color','k','position',[time(end)+300 0])
set(gca, 'Box', 'off','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1]);           % 坐标轴颜色

yyaxis left
hold on
f3 = fill([time(1:240);flipud(time(1:240))],[MA_m'+2*MA_std';flipud(MA_m'-2*MA_std')],[.7 .7 .7])
l3 = plot(time(1:240),MA_m,'color',[2 38 62]/256,'linestyle','-','linewidth',1.3)
set(f3,'EdgeColor','none','FaceAlpha',0.5)
B=datestr(time(6+12:48:end),'yyyy');
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end-2)],...
         'ylim',[-6e6 6e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-6e6:6e6:6e6,'yticklabel',-6:6:6);
y1 = ylabel({'Interdecadal component';'[10^6 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal','color','k','position',[time(1)-300 0])
% xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
% set(LG1,'location','northeast','NumColumns',5)
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
set(s1,'pos',[0.11 0.572 0.81 0.388])
set(gca,'layer','top')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t1 = text(time(130),2.5e6,['R = ',sprintf('%5.2f',r_MA(2,1)),', {\itP} < 0.01 '],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*1.2e7-0.6e7])

ax3 = axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s2=subplot(2,1,2)
hold on

f4 = fill([time(1:240);flipud(time(1:240))],[Tt_m+2*Tt_std';flipud(Tt_m-2*Tt_std')],[.7 .7 .7])
l4 = plot(time(1:240),Tt_m,'color',[2 38 62]/256,'linestyle','-','linewidth',1.3)
set(f4,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...             % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end-2)],...
         'ylim',[2.0e7 2.5e7],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',2.0e7:0.5e7:2.5e7,'yticklabel',2.0:0.5:2.5);
ylabel({'Residual component';'[10^7 km^2]'},'FontName','Arial','fontsize',12,'fontweight','normal','position',[time(1)-300 2.25e7])
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*0.5e7+2e7])
set(s2,'pos',[0.11 0.11 0.81 0.388])
% t3 = text(time(174),0.3,[sprintf('%6.4f',Tt_c_yr),' (m/s)/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time(110),2.1e7,['4.16(±1.86) × 10^4 km^2/yr, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

print(figure(11),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_6_1'],'-dpng','-r1200')
%%

figure(12)
set(gcf,'pos',[250 200 750 680])
set(gcf,'color',[1 1 1])

s1 = subplot(3,1,1)
hold on

plot(1:1440,S_ext,'color',C2,'linestyle','-','linewidth',1.3)
plot(601:840,S_crt(1:240),'color',[2 38 62]/256,'linestyle','-','linewidth',1.3)
set(gca, 'Box', 'off','layer','top','color','none', ...             % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[1 1440],...
         'ylim',[0 4e7],...
         'Xtick',0:120:1440,'xticklabel',0:120:1440,...
         'Ytick',0:2e7:4e7,'yticklabel',0:2:4);
ylabel({'Extended data';'[10^7 km^2]'},'FontName','Arial','fontsize',12,'fontweight','normal','position',[-52 2e7])
xlabel('X(t)','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[52 0.85*4e7])
set(gca,'pos',[0.12 0.72 0.79 0.25])
% t3 = text(time(174),0.3,[sprintf('%6.4f',Tt_c_yr),' (m/s)/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
LG1 = legend('Extended data','Original data');
set(LG1,'Location','southeast')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s2=subplot(3,1,2)
hold on
yyaxis right
a1 = area(time_pdo,pdo_p,'LineWidth',1,'FaceColor',C1,'EdgeColor','none','FaceAlpha',.6,'EdgeAlpha',1,'ShowBaseLine','off');
hold on
a2 = area(time_pdo,pdo_n,'LineWidth',1,'FaceColor',C2,'EdgeColor','none','FaceAlpha',.6,'EdgeAlpha',1,'ShowBaseLine','off');
xlim([time(1) time(end-2)])
ylim([-3 3])

set(gca,'Xtick',time(18:48:end),'xticklabel',[],'fontsize',12,'fontweight','normal');
set(gca,'Ytick',-3:3:3,'yticklabel',-3:3:3,'fontsize',12,'fontweight','normal');
ylabel('PDO index','Rotation',90,'FontName','Arial','fontsize',11,...
    'fontweight','normal','color','k','position',[time(end)+300 0])
set(gca, 'Box', 'off','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1]);           % 坐标轴颜色

yyaxis left
hold on
f3 = fill([time(1:240);flipud(time(1:240))],[MA_m'+2*MA_std';flipud(MA_m'-2*MA_std')],[.7 .7 .7])
l3 = plot(time(1:240),MA_m,'color',[2 38 62]/256,'linestyle','-','linewidth',1.3)
set(f3,'EdgeColor','none','FaceAlpha',0.5)
B=datestr(time(6+12:48:end),'yyyy');
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end-2)],...
         'ylim',[-6e6 6e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-6e6:6e6:6e6,'yticklabel',-6:6:6);
y1 = ylabel({'Interdecadal component';'[10^6 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal','color','k','position',[time(1)-300 0])
% xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
% set(LG1,'location','northeast','NumColumns',5)
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
set(s2,'pos',[0.12 0.385 0.79 0.25])
set(gca,'layer','top')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
t1 = text(time(130),2.5e6,['R = ',sprintf('%5.2f',r_MA(2,1)),', {\itP} < 0.01 '],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(end)-370 0.85*1.2e7-0.6e7])

ax2 = axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s3=subplot(3,1,3)
hold on

f4 = fill([time(1:240);flipud(time(1:240))],[Tt_m+2*Tt_std';flipud(Tt_m-2*Tt_std')],[.7 .7 .7])
l4 = plot(time(1:240),Tt_m,'color',[2 38 62]/256,'linestyle','-','linewidth',1.3)
set(f4,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...             % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end-2)],...
         'ylim',[2.0e7 2.5e7],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',2.0e7:0.5e7:2.5e7,'yticklabel',2.0:0.5:2.5);
ylabel({'Residual component';'[10^7 km^2]'},'FontName','Arial','fontsize',12,'fontweight','normal','position',[time(1)-300 2.25e7])
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*0.5e7+2e7])
set(s3,'pos',[0.12 0.085 0.79 0.25])
% t3 = text(time(174),0.3,[sprintf('%6.4f',Tt_c_yr),' (m/s)/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');
t4 = text(time(110),2.1e7,['4.16(±1.86) × 10^4 km^2/yr, {\itP} < 0.05'],...
    'color','r','FontName','Arial','fontsize',15,'fontweight','normal');
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

print(figure(12),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_6_2'],'-dpng','-r1200')



