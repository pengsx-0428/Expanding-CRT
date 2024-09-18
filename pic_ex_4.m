
clear;clc

load('E:\进展\文章相关\EEP_HCT\nc_reviewer\code\S_crt.mat');
S_crt = S_crt(1:240);
time = time(1:240);
Nstd = 0.01;   %Nstd:0.01-0.4(0.01,0.05,0.1,0.2);高斯白噪声标准差
NE = 800;       %NE：50/100；添加噪声的次数
Y=S_crt;

xsize=length(Y);
dd=1:1:xsize;
Ystd=std(Y);Y=Y/Ystd;
TNM=fix(log2(xsize))-1;
TNM2=TNM+2;
for kk=1:1:TNM2,
    for ii=1:1:xsize,
        allmode(ii,kk)=0.0;
    end
end
for iii=1:1:NE,
    for i=1:xsize,
        temp=randn(1,1)*Nstd;
        X1(i)=Y(i)+temp;
    end
    for jj=1:1:xsize,
        mode(jj,1) = Y(jj);
    end
    xorigin = X1;
    xend = xorigin;
    nmode = 1;
    while nmode <= TNM,
        xstart = xend;
        iter = 1;
        while iter<=10,
            [spmax, spmin, flag]=extrema(xstart);
            upper= spline(spmax(:,1),spmax(:,2),dd);
            lower= spline(spmin(:,1),spmin(:,2),dd);
            mean_ul = (upper + lower)/2;
            xstart = xstart - mean_ul;
            iter = iter +1;
        end
        xend = xend - xstart;
        nmode=nmode+1;
        for jj=1:1:xsize,
            mode(jj,nmode) = xstart(jj);
        end
    end
    for jj=1:1:xsize,
        mode(jj,nmode+1)=xend(jj);
    end
    allmode(:,:,iii)=mode;
end
% allmode=allmode/NE;
allmode=allmode*Ystd;

St_ee_all = squeeze(allmode(:,2,:) + allmode(:,3,:));
It_ee_all = squeeze(allmode(:,4,:) + allmode(:,5,:) + allmode(:,6,:) + allmode(:,7,:));
Tt_ee_all = squeeze(allmode(:,8,:));
St_ee = nanmean(St_ee_all,2);
It_ee = nanmean(It_ee_all,2);
Tt_ee = nanmean(Tt_ee_all,2);
St_ee_std = std(St_ee_all')';
It_ee_std = std(It_ee_all')';
Tt_ee_std = std(Tt_ee_all')';

Tt_m_2 = Tt_ee - detrend(Tt_ee);
Tt_yr_2 = (Tt_m_2(end) - Tt_m_2(1))/240*12
for i=1:800
    Tt_m_s = Tt_ee_all(:,i) - detrend(Tt_ee_all(:,i));
    Tt_m_all(i) = (Tt_m_s(end) - Tt_m_s(1))/240*12;
end
Tt_s_std = std(Tt_m_all)


figure(1)
set(gcf,'pos',[2650 200 750 450])
set(gcf,'color',[1 1 1])

hold on
f4 = fill([time;flipud(time)],[Tt_ee+Tt_ee_std;flipud(Tt_ee-Tt_ee_std)],[.7 .7 .7]);
l4 = plot(time,Tt_ee,'color','r','linewidth',1.3);
set(f4,'EdgeColor','none','FaceAlpha',0.5);
% l44 = plot(time,mean(Tt,2),'color',[2 38 62]/256,'linewidth',1.3)





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



%%
figure(1)
set(gcf,'pos',[2650 200 750 550])
set(gcf,'color',[1 1 1])

s2=subplot(3,1,1)
hold on
f2 = fill([time;flipud(time)],[St_ee+2*St_ee_std;flipud(St_ee-2*St_ee_std)],[.7 .7 .7])
l2 = plot(time,St_ee,'color','r','linewidth',1.3)
set(f2,'EdgeColor','none','FaceAlpha',0.5)
l22 = plot(time,mean(St,2),'color',[2 38 62]/256,'linewidth',1.3)
B=datestr(time(6+12:48:end),'yyyy');
set(gca, 'Box', 'off','layer','top','color','none', ...        % 边框
         'YAxisLocation','left',...
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-14e6 14e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-1.4e7:1.4e7:1.4e7,'yticklabel',-1.4:1.4:1.4);
ylabel({'Seasonal component';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal','position',[time(1)-500 0]);
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% text(time_chl(9),-5e6+5/6*1e7,'b','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'pos',[0.12 0.72 0.75 0.245])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(a)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*2.8e7-1.4e7])
ax2 = axes('Position',get(s2,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax2,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
% set(ax1,'pos',[0.075 0.38 0.83 0.256]);


s3=subplot(3,1,2)
hold on

f3 = fill([time;flipud(time)],[It_ee+2*It_ee_std;flipud(It_ee-2*It_ee_std)],[.7 .7 .7])
l3 = plot(time,It_ee,'color','r','linestyle','-','linewidth',1.3)
set(f3,'EdgeColor','none','FaceAlpha',0.5)
l33 = plot(time,mean(It,2),'color',[2 38 62]/256,'linewidth',1.3)
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'YAxisLocation','right',...
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],'YColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[-14e6 14e6],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',-1.4e7:1.4e7:1.4e7,'yticklabel',-1.4:1.4:1.4);
y1 = ylabel({'Interannual component';'[10^7 km^2]'},'Rotation',90,'FontName','Arial',...
    'fontsize',11,'fontweight','normal','color','k','position',[time(end)+500 0])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% text(time_chl(9),-12e6+5/6*2.4e7,'c','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'pos',[0.12 0.41 0.75 0.245])
set(gca,'layer','top')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(b)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(end)-295 0.85*2.8e7-1.4e7])
ax3 = axes('Position',get(gca,'Position'),'XAxisLocation','top','YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
set(ax3,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
% set(ax1,'pos',[0.075 0.38 0.83 0.256]);

s4=subplot(3,1,3)
hold on
f4 = fill([time;flipud(time)],[Tt_ee+2*Tt_ee_std;flipud(Tt_ee-2*Tt_ee_std)],[.7 .7 .7])
l4 = plot(time,Tt_ee,'color','r','linewidth',1.3)
set(f4,'EdgeColor','none','FaceAlpha',0.5)
l44 = plot(time,mean(Tt,2),'color',[2 38 62]/256,'linewidth',1.3)
set(gca, 'Box', 'off','layer','top','color','none', ...                                         % 边框
         'linewidth', 1.1,...
         'XGrid', 'off', 'YGrid', 'off', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.006 .006], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'XColor', [.1 .1 .1],...                           % 坐标轴颜色
         'xlim',[time(1) time(end)],...
         'ylim',[1.4e7 3.2e7],...
         'Xtick',time(6+12:48:end),'xticklabel',B,...
         'Ytick',1.4e7:9e6:3.2e7,'yticklabel',1.4:0.9:3.2);
y1=ylabel({'Residual component';'[10^7 km^2]'},'FontName','Arial',...
    'fontsize',11,'fontweight','normal','position',[time(1)-500 2.3e7])
xlabel('Year','FontName','Arial','fontsize',12,'fontweight','normal')
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
% text(time_chl(9),2e7+5/6*6e6,'d','FontName','Arial','fontsize',14,'fontweight','bold')
set(gca,'pos',[0.12 0.1 0.75 0.245])
set(gca,'FontName','Arial','fontsize',12,'fontweight','normal')
title('(c)','FontName','Arial','fontsize',13,'fontweight','bold','position',[time(1)+265 0.85*1.8e7+1.4e7])
% t1 = text(time(150),2.15e7,[sprintf('%4.2f',1.31),'×10^5 km^2/yr'],'color','r',...
%     'FontName','Arial','fontsize',15,'fontweight','normal');
t2 = text(time(90),1.7e7,['1.31(±2.46) × 10^5 km^2/yr, {\itP} < 0.05'],'color','r',...
    'FontName','Arial','fontsize',15,'fontweight','normal');
ax4 = axes('Position',get(s4,'Position'),'XAxisLocation','top','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k');
set(ax4,'XTick', [],'YTick', [],'linewidth', 1.1,'layer','top');   % 去掉xy轴刻度
% set(ax1,'pos',[0.075 0.38 0.83 0.256]);


print(figure(1),['E:\进展\文章相关\EEP_HCT\pic\v4\','pic_ex_4_4'],'-dpng','-r1200')

