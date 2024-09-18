
clear;clc;

load('E:\进展\文章相关\EEP_HCT\nc_reviewer\code\S_crt.mat');
S_crt = S_crt(1:240);
time = time(1:240);

filename =['E:\data\ENSO\MEI.txt']; %grd 文件名

%%%generate casename_grd.dat
fid=fopen(filename,'r');
% temp=fgetl(fid);
temp = fgetl(fid);
loc_time=strread(temp,'%f','delimiter',' ');
n_time=loc_time(2)-loc_time(1)+1;
for i=1:n_time
    str2=fgetl(fid);
    A(i,:)=strread(str2,'%f','delimiter',' ');
end
fclose(fid);

for i=1:n_time
    for j=1:12
        time_mei((i-1)*12+j)=datenum([sprintf('%4d',A(i,1)),'-',sprintf('%2d',j)],'yyyy-mm');
    end
end

MEI = reshape(A(:,2:end)',[44*12,1]);
MEI(MEI<-50) = nan;
C1 = [255, 128, 128]/256;
C2 = [128, 212, 255]/256;
MEI_p = MEI;MEI_p(MEI<0)=0;
MEI_n = MEI;MEI_n(MEI>0)=0;
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

Tt_yr = (Tt_m(end) - Tt_m(1))/240*12
Tt_m_2 = Tt_m - detrend(Tt_m);
Tt_yr_2 = (Tt_m_2(end) - Tt_m_2(1))/240*12

% p = polyfit(time,Tt_m',1);
% 
% Lxx = sum((time - mean(time)).^2);
% y_hat = p(1)*time + p(2);
% RSS = sum((Tt_m' - y_hat).^2);
% S = sqrt(RSS/(240-2));
% SEB = S/sqrt(Lxx);
% t_stat = p(1)/SEB;
% p_value = 2 * (1-tcdf(abs(t_stat),240-2));
p = polyfit(time,S_crt,1);

Lxx = sum((time - mean(time)).^2);
y_hat = p(1)*time + p(2);
RSS = sum((S_crt - y_hat).^2);
S = sqrt(RSS/(240-2));
SEB = S/sqrt(Lxx);
t_stat = p(1)/SEB;
p_value = 2 * (1-tcdf(abs(t_stat),240-2));

for i=1:5
    Tt_m_s = Tt(:,i) - detrend(Tt(:,i));
    Tt_m_all(i) = (Tt_m_s(end) - Tt_m_s(1))/240*12;
end
Tt_s_std = std(Tt_m_all);

[p,s] = polyfit(time,S_crt,1);
[y_fit,delta] = polyval(p,time,s);

[~,month,~] = datevec(time);

for i =1:12
    loc = find(month==i);
    cS(i) = mean(St_m(loc));
    cS_std(i) = std(St_m(loc));
end
[r,p] = corrcoef(It_m,MEI(284:523))
%%
figure(12)
set(gcf,'pos',[2650 200 750 800])
set(gcf,'color',[1 1 1])

s1=subplot(4,1,1)
hold on
l1 = plot(time,S_crt,'color',[2 38 62]/256,'linewidth',1.3)
plot(time,y_fit,'color','r','linewidth',1.3)
% plot(time,y_fit+2*delta,'r--',time,y_fit-2*delta,'r--')
B=datestr(time(6+12:48:end),'yyyy');
D = [' 0   ';'  2  ';' 4   '];
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[0 4e7],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',0:2e7:4e7,'yticklabel',D);
ylabel({'CRT area (original data)';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal','position',[time(1)-500 2e7])
set(gca,'pos',[0.12 0.77 0.75 0.19])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*4e+07])
t1 = text(time(45),1.3e7,['-2.25 × 10^4 km^2/yr, {\itP} = 0.58'],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');

ax1 = axes('Position',get(s1,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax1,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
% set(ax1,'pos',[0.075 0.38 0.83 0.256]);

s2=subplot(4,1,2)
hold on
f2 = fill([time;flipud(time)],[St_m'+2*St_std';flipud(St_m'-2*St_std')],[.7 .7 .7])
l2 = plot(time,St_m,'color',[2 38 62]/256,'linewidth',1.3)
set(f2,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...        % 边框
         'YAxisLocation','right',...
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-14e6 14e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-1.4e7:1.4e7:1.4e7,'yticklabel',-1.4:1.4:1.4);
ylabel({'Seasonal component';'[10^7 km^2]'},'Rotation',90,'FontName','Arial',...
    'fontsize',11,'fontweight','normal','position',[time(end)+500 4]);
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% text(time_chl(9),-5e6+5/6*1e7,'b','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'pos',[0.12 0.54 0.75 0.19])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(end)-295 0.85*2.8e7-1.4e7])
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
% set(ax1,'pos',[0.075 0.38 0.83 0.256]);


s3=subplot(4,1,3)
hold on
yyaxis right
area(time_mei,MEI_p,'LineWidth',1,'FaceColor',C1,'EdgeColor','none','FaceAlpha',.8,'EdgeAlpha',1,'ShowBaseLine','off');
hold on
area(time_mei,MEI_n,'LineWidth',1,'FaceColor',C2,'EdgeColor','none','FaceAlpha',.8,'EdgeAlpha',1,'ShowBaseLine','off');
xlim([time(1) time(end)])
ylim([-3 3])
set(gca,'Xtick',time(18:48:end),'xticklabel',[],'fontsize',12,'fontweight','normal');
set(gca,'Ytick',-3:3:3,'yticklabel',-3:3:3,'fontsize',12,'fontweight','normal');
ylabel('MEI','Rotation',90,'FontName','Arial','fontsize',11,...
    'fontweight','normal','color','k','position',[time(end)+500 0])
set(gca, 'Box', 'off','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1]);           % 坐标轴颜色

yyaxis left
f3 = fill([time;flipud(time)],[It_m'+2*It_std';flipud(It_m'-2*It_std')],[.7 .7 .7])
l3 = plot(time,It_m,'color',[2 38 62]/256,'linestyle','-','linewidth',1.3)
set(f3,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-14e6 14e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-1.4e7:1.4e7:1.4e7,'yticklabel',-1.4:1.4:1.4);
y1 = ylabel({'Interannual component';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal','color','k','position',[time(1)-500 11.444])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% text(time_chl(9),-12e6+5/6*2.4e7,'c','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'pos',[0.12 0.31 0.75 0.19])
set(gca,'layer','top')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*2.8e7-1.4e7])
ax3 = axes('Position',get(s3,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度

s4=subplot(4,1,4)
hold on
f4 = fill([time;flipud(time)],[Tt_m'+2*Tt_std';flipud(Tt_m'-2*Tt_std')],[.7 .7 .7])
l4 = plot(time,Tt_m,'color',[2 38 62]/256,'linewidth',1.3)
set(f4,'EdgeColor','none','FaceAlpha',0.5)
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'YAxisLocation','right',...
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[2e7 3e7],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',2e7:5e6:3e7,'yticklabel',2:0.5:3);
y1=ylabel({'Residual component';'[10^7 km^2]'},'Rotation',90,'FontName','Arial',...
    'fontsize',11,'fontweight','normal','position',[time(end)+500 2.5e7])
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% text(time_chl(9),2e7+5/6*6e6,'d','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'pos',[0.12 0.08 0.75 0.19])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(d)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(end)-295 0.85*1e7+2e7])
% t1 = text(time_chl(150),2.2e7,[sprintf('%4.2f',8.09),'×10^4 km^2/yr'],'color','k',...
%     'FontName','Arial','fontsize',12,'fontweight','normal');

t2 = text(time(30),2.7e7,['1.87(±0.82) × 10^5 km^2/yr, {\itP} < 0.05'],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
ax4 = axes('Position',get(s4,'Position'),'XAxisLocation','top','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax4,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
% set(ax1,'pos',[0.075 0.38 0.83 0.256]);

% set(gca,'pos',[0.12 0.54 0.75 0.19])
ax5 = axes('pos',[0.6 0.55 0.3 0.07],'XAxisLocation','bottom','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax5,'pos',[0.5 0.535 0.33 0.095])
set(ax5,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
set(ax5,'XColor','none')
set(ax5,'Ycolor','none')
hold on
E1 = errorbar(1:12,cS,cS_std);

set(E1,  'LineStyle', '-', 'Color', 'k',...
         'LineWidth', 1.1, 'Marker', 'o', 'MarkerSize', 3, ...
         'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , 'k')
BB = datestr(time(6:17),'mmm');
D = [' 1   ';'  2  ';' 3   '];
set(gca, 'Box', 'off','layer','top','color','none', ...              % 边框
         'linewidth', 1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[-1.5 14.5],...
         'ylim',[-6e6 6e6],...
         'Xtick',1:12,'xticklabel',[],...
         'Ytick',-4e6:4e6:4e6,'yticklabel',[]);
text(-0.4,-1.5e6,'Jan','FontName','Arial','fontsize',14,'fontweight','normal')
text(11.6,-3.5e6,'Dec','FontName','Arial','fontsize',14,'fontweight','normal')
line([-0.6 14 14 -0.6 -0.6],[-4.9e6 -4.9e6 3.9e6 3.9e6 -4.9e6],'color','k','linewidth',1)

print(figure(12),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_2_3'],'-dpng','-r1200')




