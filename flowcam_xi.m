clear all; close all; clc; load flowcam.mat; 

%% calculate exponents & uncertainties

l = who;
for i = 1:length(l);
    a(i) = plfit(eval(l{i}),'xmin',2); % power-law exponent for each sample
    [au(i),~,~] = plvar(eval(l{i}),'xmin',2,'silent'); % bootstrap uncertainty for each sample
    i
end

%% calculate means & uncertainties

A0 = (sum(a(1:3)./au(1:3).^2))./sum(au(1:3).^-2); % (uncertainty-weighted) mean exponent -- control 
U0 = sqrt(sum((a(1:3)-A0).^2./au(1:3).^2)./sum(au(1:3).^-2))./sqrt(3); % uncertainty in the mean
A0 = A0-a(4); % mean change in slope in controls
U0 = sqrt(U0^2+au(4)^2); % uncertainty in mean change in slope in controls

A1 = (sum(a(5:7)./au(5:7).^2))./sum(au(5:7).^-2); % low treatment
U1 = sqrt(sum((a(5:7)-A1).^2./au(5:7).^2)./sum(au(5:7).^-2))./sqrt(3);
A1 = A1-a(8);
U1 = sqrt(U1^2+au(8)^2);

A2 = (sum(a(1:3)./au(1:3).^2))./sum(au(1:3).^-2); % medium treatment
U2 = sqrt(sum((a(9:11)-A2).^2./au(9:11).^2)./sum(au(9:11).^-2))./sqrt(3);
A2 = A2-a(12);
U2 = sqrt(U2^2+au(12)^2);

a3(1:3) = a(13:15)-a(16); a3(4:6) = a(21:23)-a(24); % high treatment, both stations -- note have to subtract station's t0's separately
au3(1:3) = sqrt(au(13:15).^2+au(16).^2); au3(4:6) = sqrt(au(21:23).^2+au(24).^2);
A3 = sum(a3./au3.^2)./sum(au3.^-2);
U3 = sqrt(sum((a3-A3).^2./au3.^2)./sum(au3.^-2))./sqrt(6);

a31 = a3(1:3); % high treatment, station 1
au31(1:3) = sqrt(au(13:15).^2+au(16).^2);
A31 = sum(a31./au31.^2)./sum(au31.^-2);
U31 = sqrt(sum((a31-A31).^2./au31.^2)./sum(au31.^-2))./sqrt(3);

a32 = a3(4:6); % high treatment, station 2
au32(1:3) = sqrt(au(21:23).^2+au(24).^2);
A32 = sum(a32./au32.^2)./sum(au32.^-2);
U32 = sqrt(sum((a32-A32).^2./au32.^2)./sum(au32.^-2))./sqrt(3);

%% plot

y(1:3) = a(1:3)-a(4); % individual bottles' values of \Delta \xi \pm uncertainty
y(4:6) = a(5:7)-a(8);
y(7:9) = a(9:11)-a(12);
y(10:12) = a(13:15)-a(16);
y(13:15) = a(21:23)-a(24);
yu(1:3) = sqrt(au(1:3).^2+au(4).^2)
yu(4:6) = sqrt(au(5:7).^2+au(8).^2)
yu(7:9) = sqrt(au(9:11).^2+au(12).^2)
yu(10:12) = sqrt(au(13:15).^2+au(16).^2);
yu(13:15) = sqrt(au(21:23).^2+au(24).^2);
x = [0 0 0 1 1 1 2 2 2 3 3 3 3.25 3.25 3.25];

l1 = errorbar(x,-y,yu,'o','MarkerSize',6,...
    'Color',[0 0 0],'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'linewidth',1.5)
hold on
l2 = errorbar([0 1 2 3],-[A0 A1 A2 A31],[U0 U1 U2 U31],'o','MarkerSize',12,...
    'Color',[.59 .07 .39],'MarkerEdgeColor',[.59 .07 .39],'MarkerFaceColor',[.59 .07 .39],'linewidth',3)
errorbar([3.25],-[A32],[U32],'s','MarkerSize',12,...
    'Color',[.59 .07 .39],'MarkerEdgeColor',[.59 .07 .39],'MarkerFaceColor',[.59 .07 .39],'linewidth',3)
errorbar([3.5],-[A3],[U3],'d','MarkerSize',12,...
    'Color',[.59 .07 .39],'MarkerEdgeColor',[.59 .07 .39],'MarkerFaceColor',[.59 .07 .39],'linewidth',3)
set(gca,'ticklabelinterpreter','latex','fontsize',18,'xtick',[0:3 3.25 3.5],'xticklabel',{'Control','Low','Medium','High','(Station 2)','(Pooled)'},'ytick',[-.7:.1:.1])
xtickangle(90)
axis([-.5 4 -0.75 0.15])
ylabel('Size distribution change: $\xi_f-\xi_i$','interpreter','latex')


% note p-value, slope, & intercept calculated separately
l5 = plot(linspace(-.25,3.75),.0494-.1504.*linspace(-.25,3.75),'color',[0.7373    0.5882    0.0275],'linewidth',3)
lgnd = legend([l1 l2 l5],'Bottles','Averages','$p \leq 0.030$','location','southwest') 
set(lgnd,'interpreter','latex')

clc;