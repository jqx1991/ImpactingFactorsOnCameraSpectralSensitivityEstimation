% plot spectral sensitivity curves
load('E:\Dropbox\Works\Matlab\Papers\ACA2016\TestingUsingMinorTrainingSamples_result.mat')
CSS_D65 = CSS{4};
CSS_D65 = CSS_D65/max(CSS_D65(:));
figure;
wl = [380:10:780]';
hold on;
R = plot(wl,CSS_D65(:,1),'color',[230 130 120]/255,'LineStyle','-','lineWidth',3);
G = plot(wl,CSS_D65(:,2),'color',[170 235 160]/255,'LineStyle','-','lineWidth',3);
B = plot(wl,CSS_D65(:,3),'color',[110 170 240]/255,'LineStyle','-','lineWidth',3);

xlim([380 780]);ylim([0 1.1]);
box on;
set(gca,'TickLabelInterpreter','LaTex','FontSize',14);
xlabel('$\textrm{Wavelength }(\textrm{nm})$','Interpreter','latex','FontSize',16);
ylabel('$\textrm{Relative camera spectral sensitivity }$','Interpreter','latex','FontSize',16);
set(gcf,'color','w','Units','inches','Position',[2 2 6.5 4.5]);

legend({'estimated $\mathbf{S}^{R}$',...
        'estimated $\mathbf{S}^{G}$',...
        'estimated $\mathbf{S}^{B}$',...
        },'Interpreter','latex','FontSize',14)
    
PHi = 0:0.001:1;
nl0 = real(PHi+nonlinearCoef0(1)).^(nonlinearCoef0(3)) + nonlinearCoef0(2);
nl = real(PHi+nonlinearCoef(1)).^(nonlinearCoef(3)) + nonlinearCoef(2);
figure;box on;hold on
l1 = line([0 1],[0 1],'LineWidth',2,'color',[200 200 200]/255)
l2 = plot(PHi,nl0,'LineWidth',2.5,'color',[119 163 188]/255);
l3 = plot(PHi,nl,'LineStyle','--','LineWidth',2.5,'color',[229 134 140]/255);
axis equal;xlim([0 1]);ylim([0 1]);
legend([l2,l3,l1],{...
        'initial fitting: $f(\Phi ) = (\Phi  + 0.6172)^{0.9871} - 0.6207$',...
        'optimized: $f(\Phi ) = (\Phi - 0.00139)^{0.9875} + 0.00176$',...
        '$f(\Phi ) = \Phi$',...
        },'Interpreter','latex','FontSize',11,'Box','off')
xlabel('$\Phi$','Interpreter','latex','FontSize',15);
ylabel('$f(\Phi)$','Interpreter','latex','FontSize',15,'Rotation',0);

set(gca,'TickLabelInterpreter','LaTex','FontSize',13,'LineWidth',1.5);
set(gca,'XTick',[0:0.2:1],'YTick',[0:0.2:1])
set(gcf,'color','w','Units','inches','Position',[3 3 6.5 5.5]);


