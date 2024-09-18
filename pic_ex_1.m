
clear;clc;

load('E:\进展\文章相关\EEP_HCT\nc_reviewer\code\S_crt.mat');
S_crt = S_crt(1:240);
time = time(1:240);

%%

[imf_2, residual_2] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',2);
[imf_3, residual_3] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',3);
[imf_4, residual_4] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',4);
[imf_5, residual_5] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',5);
[imf_6, residual_6] = emd(S_crt,'SiftRelativeTolerance',0,'SiftMaxIterations',6);

imf1 = [imf_2(:,1),imf_3(:,1),imf_4(:,1),imf_5(:,1),imf_6(:,1)];
imf2 = [imf_2(:,2),imf_3(:,2),imf_4(:,2),imf_5(:,2),imf_6(:,2)];
imf3 = [imf_2(:,3),imf_3(:,3),imf_4(:,3),imf_5(:,3),imf_6(:,3)];
imf4 = [imf_2(:,4),imf_3(:,4),imf_4(:,4),imf_5(:,4),imf_6(:,4)];
imf5 = [imf_2(:,5),imf_3(:,5),imf_4(:,5),imf_5(:,5),imf_6(:,5)];
imf6 = [imf_2(:,6),imf_6(:,6)];

imf1_m = mean(imf1,2);
imf2_m = mean(imf2,2);
imf3_m = mean(imf3,2);
imf4_m = mean(imf4,2);
imf5_m = mean(imf5,2);
imf6_m = mean(imf6,2);

imf1_std = std(imf1');
imf2_std = std(imf2');
imf3_std = std(imf3');
imf4_std = std(imf4');
imf5_std = std(imf5');
imf6_std = std(imf6');

res = [residual_2,residual_3,residual_4,residual_5,residual_6];
res_m = mean(res,2);
res_std = std(res');

%%

figure(1)
set(gcf,'pos',[2650 200 750 800])
set(gcf,'color',[1 1 1])

% s1=subplot(8,1,1)
% l1 = plot(time,S_crt,'color',[2 38 62]/256,'linewidth',1.2)
% B=datestr(time(6+12:48:end),'yyyy');
% D = [' 0   ';'  2  ';' 4   '];
% set(gca, 'Box', 'on','layer','top', ...                 % 边框
%          'linewidth', 1.1,...
%          'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
%          'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
%          'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
%          'XColor', [.1 .1 .1],...                           % 坐标轴颜色
%          'xlim',[time(1) time(end)],...
%          'ylim',[0 4e7],...
%          'Xtick',time(6+12:48:end),'xticklabel',B,...
%          'Ytick',0:2e7:4e7,'yticklabel',D);
% ylabel({'Signal';'[10^7 km^2]'},'FontName','Arial',...
%     'fontsize',11,'fontweight','normal')
% set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(1)+265 0.85*4e+07])

s2=subplot(7,1,1)
hold on
f2 = fill([time;flipud(time)],[imf1_m+2*imf1_std';flipud(imf1_m-2*imf1_std')],[.7 .7 .7])
l2 = plot(time,imf1_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f2,'EdgeColor','none','FaceAlpha',0.5)
B=datestr(time(6+12:48:end),'yyyy');
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-10e6 10e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-10e6:10e6:10e6,'yticklabel',-1:1:1);
ylabel({'IMF1';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

s3=subplot(7,1,2)
hold on
f3 = fill([time;flipud(time)],[imf2_m+2*imf2_std';flipud(imf2_m-2*imf2_std')],[.7 .7 .7])
l3 = plot(time,imf2_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f3,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-10e6 10e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-10e6:10e6:10e6,'yticklabel',-1:1:1);
ylabel({'IMF2';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

s4=subplot(7,1,3)
hold on
f4 = fill([time;flipud(time)],[imf3_m+2*imf3_std';flipud(imf3_m-2*imf3_std')],[.7 .7 .7])
l4 = plot(time,imf3_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f4,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-10e6 10e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-10e6:10e6:10e6,'yticklabel',-1:1:1);
ylabel({'IMF3';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

s5=subplot(7,1,4)
hold on
f5 = fill([time;flipud(time)],[imf4_m+2*imf4_std';flipud(imf4_m-2*imf4_std')],[.7 .7 .7])
l5 = plot(time,imf4_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f5,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-10e6 10e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-10e6:10e6:10e6,'yticklabel',-1:1:1);
ylabel({'IMF4';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

s6=subplot(7,1,5)
hold on
f6 = fill([time;flipud(time)],[imf5_m+2*imf5_std';flipud(imf5_m-2*imf5_std')],[.7 .7 .7])
l6 = plot(time,imf5_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f6,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-10e6 10e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-10e6:10e6:10e6,'yticklabel',-1:1:1);
ylabel({'IMF5';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

s7=subplot(7,1,6)
hold on
f7 = fill([time;flipud(time)],[imf6_m+2*imf6_std';flipud(imf6_m-2*imf6_std')],[.7 .7 .7])
l7 = plot(time,imf6_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f7,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-10e5 10e5],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-10e5:10e5:10e5,'yticklabel',-1:1:1);
ylabel({'IMF6';'[10^6 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

s8=subplot(7,1,7)
hold on
f8 = fill([time;flipud(time)],[res_m+2*res_std';flipud(res_m-2*res_std')],[.7 .7 .7])
l8 = plot(time,res_m,'color',[2 38 62]/256,'linewidth',1.2)
set(f8,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'on','layer','top', ...        % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'in', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[2e7 3e7],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',2e7:5e6:3e7,'yticklabel',2:0.5:3);
ylabel({'Residual';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal');
xlabel('Year','FontName','Arial','fontsize',11,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',11,'fontweight','normal')
% title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time_chl(end)-295 0.85*1e7-0.5e7])

% set(s1,'pos',[0.12 0.89 0.85 0.09])
set(s2,'pos',[0.12 0.06+0.82/6*6 0.85 0.1])
set(s3,'pos',[0.12 0.06+0.82/6*5 0.85 0.1])
set(s4,'pos',[0.12 0.06+0.82/6*4 0.85 0.1])
set(s5,'pos',[0.12 0.06+0.82/6*3 0.85 0.1])
set(s6,'pos',[0.12 0.06+0.82/6*2 0.85 0.1])
set(s7,'pos',[0.12 0.06+0.82/6*1 0.85 0.1])
set(s8,'pos',[0.12 0.06 0.85 0.1])

print(figure(1),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_1'],'-dpng','-r1200')

