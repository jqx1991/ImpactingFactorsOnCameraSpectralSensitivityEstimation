load('E:\Dropbox\Works\Matlab\Papers\ACA2016\TestingUsingMinorTrainingSamples_A_result.mat');
for i = 1:10
    DE_D65(i,:) = mean(DeltaE_final_D65{i});
    err_D65(i) = std(DeltaE_final_D65{i})/sqrt(length(DeltaE_final_D65{i}));
    ts_D65(i,:) = tinv([0.05 0.95],length(DeltaE_final_D65{i})-1);
    CI_D65(i,:) = ts_D65(i,:)*err_D65(i);
    
    DE_A(i) = mean(DeltaE_final_A{i});
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
num = 4:2:22;
hold on;
for i = 1:10
    bar(num(i),DE_A(i),1.5,'FaceColor', cmap(i,:));
end
errorbar(num, DE_A, CI_A(:,1),CI_A(:,2), 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k');
set(gcf,'color','w','Units','inches','Position',[2 2 6.5 4.5]);
set(gca,'TickLabelInterpreter','LaTex','FontSize',12);
set(gca,'xtick',[4:2:22],'ytick',[0:0.5:6]);
xlim([2 24]);ylim([0 6])
xlabel('$\textrm{The number of training samples}$','Interpreter','latex','FontSize',14);
ylabel('$\Delta E_{00} \textrm{ of testing samples}$','Interpreter','latex','FontSize',14);
box on
