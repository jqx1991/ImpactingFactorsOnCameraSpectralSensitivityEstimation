% Using CSS obtained under D65 to test illuminant A
load('E:\Dropbox\Works\Matlab\Papers\ACA2016\TestingUsingMinorTrainingSamples_result.mat')

load('M:\D3x\A\data\RGB_mean_ranked.mat')
RGB_A = RGB_mean_ranked;
RGB_A = RGB_A(1:24,:);
load('E:\Dropbox\Works\papers\ResponsePrediction\SPD_Central_A.mat')
for i = 1:24
    SPD_A(i,:) = (SPD_Central_A(2*i-1,:)+SPD_Central_A(2*i,:))/2;
end
SPD_A = SPD_A(:,1:10:end);
clear RGB_mean_ranked SPD_Central_A

const_geometry = (pi/4)*((1/4)^2); % #F = 4
ISO_A = 1*ones(24,1);
ExposureTime_A = 1/15*ones(24,1);
Const_A = ISO_A.*ExposureTime_A * const_geometry * 10;
RGB_reconst_test_A = (real( ( diag(Const_A)*SPD_A*CSS{3}*CrossTalkMtx{3} + nonlinearCoef{3}(1) ).^nonlinearCoef{3}(3) )) + nonlinearCoef{3}(2);
DeltaE_final_A = sRGB2CIEDeltaE(RGB_A.^(1/2.2),RGB_reconst_test_A.^(1/2.2),'cie00');

load('E:\Dropbox\Works\Matlab\Ruixinwei\ToolFunctions\SpectralData\ClassicColorCheckerSpectralReflectance.mat')
sRGB = SpecRef2sRGB(ClassicColorCheckerSpectralReflectance',380,780);
sRGB(sRGB>1) = 1;
% training samples
figure;
for i = 1:24
    hold on;
    bar(i,DeltaE_final_A(i),1,'facecolor',sRGB(i,:));
end
hold on;
l = line([0 50],[mean(DeltaE_final_A) mean(DeltaE_final_A)],'LineStyle','--','LineWidth',2,'Color',[160 160 160]/255);
xlim([0 25]);ylim([0 3]);
set(gca,'TickLabelInterpreter','LaTex','FontSize',12);
set(gca,'xtick',[1,5:5:20,24],'ytick',[0:0.5:3])
xlabel('$\textrm{Color patch No. }$','Interpreter','latex','FontSize',14);
ylabel('$\Delta E_{00} \textrm{ of testing samples}$','FontSize',14);
set(gcf,'color','w','Units','inches','Position',[2 2 6.5 4.5]);
box on;
legend(l,{['$\textrm{mean }\Delta {E_{00}}: $ ',num2str(mean(DeltaE_final_A))]},'Interpreter','latex','FontSize',12,'Box','off')



% Using CSS obtained under illuminant A to test D65
load('E:\Dropbox\Works\Matlab\Papers\ACA2016\TestingUsingMinorTrainingSamples_A_result.mat')

load('M:\D3x\Central\data\RGB_mean_ranked.mat')
RGB_D65 = RGB_mean_ranked;
RGB_D65 = RGB_D65(1:96,:);
RGB_D65 = RGB_D65([4:9,16:21,28:33,40:45],:);
load('E:\Dropbox\Works\papers\ResponsePrediction\SPD_Central.mat')
for i = 1:96
    SPD_D65(i,:) = (SPD_Central(2*i-1,:)+SPD_Central(2*i,:))/2;
end
SPD_D65 = SPD_D65([4:9,16:21,28:33,40:45],1:10:end);
clear RGB_mean_ranked SPD_Central

const_geometry = (pi/4)*((1/4)^2); % #F = 4
ISO_D65 = 1*ones(24,1);
ExposureTime_D65 = 1/15*ones(24,1);
Const_D65 = ISO_D65.*ExposureTime_D65 * const_geometry * 10;
RGB_reconst_test_D65 = (real( ( diag(Const_D65)*SPD_D65*CSS{5}*CrossTalkMtx{5} + nonlinearCoef{5}(1) ).^nonlinearCoef{5}(3) )) + nonlinearCoef{5}(2);
DeltaE_final_D65 = sRGB2CIEDeltaE(RGB_D65.^(1/2.2),RGB_reconst_test_D65.^(1/2.2),'cie00');

load('E:\Dropbox\Works\Matlab\Ruixinwei\ToolFunctions\SpectralData\ClassicColorCheckerSpectralReflectance.mat')
sRGB = SpecRef2sRGB(ClassicColorCheckerSpectralReflectance',380,780);
sRGB(sRGB>1) = 1;
figure;
for i = 1:24
    hold on;
    bar(i,DeltaE_final_D65(i),1,'facecolor',sRGB(i,:));
end
hold on;
l = line([0 50],[mean(DeltaE_final_D65) mean(DeltaE_final_D65)],'LineStyle','--','LineWidth',2,'Color',[160 160 160]/255);
xlim([0 25]);ylim([0 3]);
set(gca,'TickLabelInterpreter','LaTex','FontSize',12);
set(gca,'xtick',[1,5:5:20,24],'ytick',[0:0.5:3])
xlabel('$\textrm{Color patch No. }$','Interpreter','latex','FontSize',14);
ylabel('$\Delta E_{00} \textrm{ of testing samples}$','FontSize',14);
set(gcf,'color','w','Units','inches','Position',[2 2 6.5 4.5]);
box on;
legend(l,{['$\textrm{mean }\Delta {E_{00}}: $ ',num2str(mean(DeltaE_final_D65))]},'Interpreter','latex','FontSize',12,'Box','off')
