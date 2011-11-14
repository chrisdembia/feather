function qqktplots( fname )

% load file

T = load([fname '_t.mat'],'-ascii');
Q = load([fname '_q.mat'],'-ascii');
QK = load([fname '_qk.mat'],'-ascii');
Tau = load([fname '_tau.mat'],'-ascii');

VP = load('path4accordion.mat','-ascii');

idx = find( T > 6,1,'first');
T = T(1:idx,1);
size(T)
Q = Q(1:idx,:);
size(Q)
QK = QK(1:idx,:);
size(QK)
Tau = Tau(1:idx,:);
size(Tau)
% VP = VP(1:idx,:);
[maxTidx N] = size(Q);
close all
% figure
% hold on

for i = 1:N
%     
%     if i <=6 
%     poss(1) = .005;
% elseif i > 6
%     poss(1) = .5;
% end
% poss(2) = (1-(mod(i-1,6))/6)-1/6;
% poss(3) = .495;
% poss(4) = 1/6;
%     set(g,'Position',[poss(1) poss(2) 1/2 1/6]);
%     g = subplot('position',poss);
    h =figure;
    set(h,'Position',[100 100 740 300],'PaperPositionMode','auto');

    
      hold on;
    [ax,h1,h2] = plotyy([T T],[Q(:,i) QK(:,i)],T,Tau(:,i));
    xlabel('t (s)')
    ylabel('q (rad)')
    set(get(ax(1),'Ylabel'),'String','q (rad)')
    set(get(ax(2),'Ylabel'),'String','\tau (Nm)')
    set(ax(1),'YTick',[-pi -pi/2 0 pi/2 pi],...
              'YTickLabel',{'','','0','',''},...
              'YLim',[-pi/.9 pi/.9],'XLim',[0 6]);
          set(ax(2),'XLim',[0 6]);

%     set(ax(2),'YTickLabel',{},'XTickLabel',{},'XLim',[0 6],'LineWidth',2,'Box','on');%,...
% %               'YLim',[-1.25*pi 1.25*pi])
%     set(ax(1), 'YLim',[-1.25*pi 1.25*pi],...
%            'XTickLabel',{},'XLim',[0 6],'LineWidth',2,'Box','on');
%        text(0,-1.25*pi,[' ' num2str(i)],'VerticalAlignment','bottom')
% hold on; axis off; box on; 
    text(0,-pi,'-\pi ','HorizontalAlignment','right');
    text(0,-pi/2,'-\pi/2 ','HorizontalAlignment','right');
    text(0,pi/2,'\pi/2 ','HorizontalAlignment','right');
    text(0,pi,'\pi ','HorizontalAlignment','right');

    legend('Q','Qk')
    
        plot(VP(:,end),VP(:,i),'r.','MarkerSize',10)
    plot([0 T(end)],[0 0],'Color',[.5 .5 .5])  
%     pause
    print(['c:\dropbox\crobobauts\proj\tex\qqkt' num2str(i) '.png'],'-dpng','-r300');
        
end
% print(['c:\dropbox\crobobauts\proj\tex\qqktall.png'],'-dpng','-r500');
% plot 11 plots
