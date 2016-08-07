load('E:\Dropbox\Works\Matlab\Papers\ACA2016\TestingUsingMinorTrainingSamples_result.mat');
for i = 1:12
    DE_D65(i,:) = mean(DeltaE_final_D65{i});
    err_D65(i) = std(DeltaE_final_D65{i})/sqrt(length(DeltaE_final_D65{i}));
    ts_D65(i,:) = tinv([0.05 0.95],length(DeltaE_final_D65{i})-1);
    CI_D65(i,:) = ts_D65(i,:)*err_D65(i);
    DE_A(i,:) = mean(DeltaE_final_A{i});
    err_A(i) = std(DeltaE_final_A{i})/sqrt(length(DeltaE_final_A{i}));
    ts_A(i,:) = tinv([0.05 0.95],length(DeltaE_final_A{i})-1);
    CI_A(i,:) = ts_A(i,:)*err_A(i);
end

I = imread('E:\index.png');
RGB = reshape(double(I(401:500,:,:)),100*1200,3);
[~,C] = kmeans(RGB,48);
clear I RGB
cmap = C;
cmap = cmap/255;

figure;
num = 4:4:48;
hold on;
for i = 1:12
    bar(num(i),DE_D65(i),3,'FaceColor', cmap(i,:));
end
errorbar(num, DE_D65, CI_D65(:,1),CI_D65(:,2), 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k');
set(gcf,'color','w','Units','inches','Position',[2 2 7 5]);
set(gca,'TickLabelInterpreter','LaTex','FontSize',12);
set(gca,'xtick',[4:4:48],'ytick',[0:0.5:6]);
xlim([0 52]);ylim([0 6])
xlabel('$\textrm{The number of training samples}$','Interpreter','latex','FontSize',14);
ylabel('$\Delta E_{00} \textrm{ of testing samples}$','Interpreter','latex','FontSize',14);
box on

load('E:\Dropbox\Works\Matlab\Papers\ACA2016\LOOCV_result.mat');
[~,order] = sort(DeltaE_test);
BestPatchOrder = order(1:24);
load('E:\Dropbox\Works\Matlab\Ruixinwei\ToolFunctions\SpectralData\DSGColorCheckerSpectralReflectance_order.mat')
sRGB = SpecRef2sRGB(DSGColorCheckerSpectralReflectance,400,780);
sRGB(sRGB>1) = 1;
PeripheralIdx = [1:14,15:14:113,28:14:126,127:140];
sRGB(PeripheralIdx,:) = [];
sRGB_best = sRGB(BestPatchOrder,:);
load('E:\Dropbox\Works\papers\ResponsePrediction\SPD_Central.mat')
for i = 1:96
    SPD_D65(i,:) = (SPD_Central(2*i-1,:)+SPD_Central(2*i,:))/2;
end
SPD_D65 = SPD_D65(:,1:5:end);
SPD_D65_best = SPD_D65(BestPatchOrder,:);
figure;
hold on;
for i = 1:24
    plot([380:5:780],SPD_D65_best(i,:),'color',sRGB_best(i,:),'LineWidth',2);
end
set(gcf,'color','w','Units','inches','Position',[2 2 6.5 4.5]);
set(gca,'TickLabelInterpreter','LaTex','FontSize',12);
set(gca,'xtick',[380:100:780],'ytick',[0:1:5]*10^-3);
xlabel('$\textrm{Wavelength }(\textrm{nm})$','Interpreter','latex','FontSize',14);
ylabel('$\textrm{The spectral radiance }(\textrm{W} \cdot \textrm{sr}^{ - 1} \cdot \textrm{m}^{-3})$','Interpreter','latex','FontSize',14);
xlim([380 780]);ylim([0 5]*10^-3);
box on
