%Mutual information quantified between the two dynamics and between the stress and Crz1 dynamics 
%MI between two dynamicsx
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_5";
folder{2}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_27_2020\GCaMP3_ScIDR_dZF_GFP_2020_09_27_3\';
folder{3}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_10_03_2020\GCaMP3_ScIDR_dZF_GFP_2020_10_03_3\';
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_5";
folder{5}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_13_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_13_3\';
folder{6}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_17_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_17_3\';
folder{7}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_5";
folder{8}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_4_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_04_3\';
folder{9}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_7_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_07_3\';
folder{10}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_5";
folder{11}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_20_2020\GCaMP3_CaIDR_dZF_GFP_2020_02_20_3\';
folder{12}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_19_2020\GCaMP3_CaIDR_dZF_GFP_2020_03_19_3\';
folder{13}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_5";
folder{14}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_10_03_2020\GCaMP3_SpIDR_dZF_GFP_2020_10_03_3\';
folder{15}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_10_11_2020\GCaMP3_SpIDR_dZF_GFP_2020_10_11_3\';
CaMI = cell(1,15);
for i = 1:15
    cd(folder{i})
    CaMI{i} = load('MI.mat');
end
%%

% MI between environments and Crz1 dynamics
cd('C:\Users\Ian\Documents\MATLAB\Crz1\MI quantification\SteadyState')
SMI = load('MIReporter.mat');

%x = [SMI.SpMI,TMI.SpMI;
%    mean([SMI.CaMI1,SMI.CaMI2,SMI.CaMI3]),mean([TMI.CaMI1,TMI.CaMI2,TMI.CaMI3]);
%    mean([SMI.KlMI1,SMI.KlMI2,SMI.KlMI3]),mean([TMI.KlMI1,TMI.KlMI2,TMI.KlMI3]);
%    mean([SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3]),mean([TMI.ZrMI1,TMI.ZrMI2,TMI.ZrMI3]);
%    SMI.ScMI,TMI.ScMI];
x = [
    mean([SMI.ScMI1,SMI.ScMI2,SMI.ScMI3]);
    mean([SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3]);
    mean([SMI.KlMI1,SMI.KlMI2,SMI.KlMI3]);
    mean([SMI.CaMI1,SMI.CaMI2,SMI.CaMI3]);
    mean([SMI.SpMI1,SMI.SpMI2,SMI.SpMI3]);];
errx = [std([SMI.ScMI1,SMI.ScMI2,SMI.ScMI3],1)...
    std([SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3],1)...
    std([SMI.KlMI1,SMI.KlMI2,SMI.KlMI3],1)...
    std([SMI.CaMI1,SMI.CaMI2,SMI.CaMI3],1)...
    std([SMI.SpMI1,SMI.SpMI2,SMI.SpMI3],1)...
    ]/3^0.5;
TTx = {[SMI.ScMI1,SMI.ScMI2,SMI.ScMI3];
    [SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3];
    [SMI.KlMI1,SMI.KlMI2,SMI.KlMI3];
    [SMI.CaMI1,SMI.CaMI2,SMI.CaMI3]; 
    [SMI.SpMI1,SMI.SpMI2,SMI.SpMI3]};
m = nan(5,1); errm = m;
%ami = mean(CaMI{(1)}.I);m(1)=mean(ami(:,1));
%ami = mean(CaMI{(11)}.I);m(5)=mean(ami(:,1));
%errm([1,5])=0;
TTm = cell(1,5);% TTm{1}=m(1);TTm{5}=m(5);
for i = 1:5
    ami1 = mean(CaMI{3*i-2}.I);
    ami2 = mean(CaMI{3*i-1}.I);
    ami3 = mean(CaMI{3*i}.I);
    m(i) = mean([mean(ami1(:,1)),mean(ami2(:,1)),mean(ami3(:,1))]);
    errm(i)=std([mean(ami1(:,1)),mean(ami2(:,1)),mean(ami3(:,1))],1)/3^0.5;
    TTm{i}=[mean(ami1(:,1)),mean(ami2(:,1)),mean(ami3(:,1))];
end

xm = [m,x];
figure(5)
bar(xm)
hold on
dif = 0.15;
er1 = errorbar([1:length(m)]'-dif,m,errm*1.96);
er1.Color = [0 0 0];
er1.LineStyle = 'none';
er1 = errorbar([1:length(x)]'+dif,x,errx*1.96);
er1.Color = [0 0 0];
er1.LineStyle = 'none';
hold off
xticklabels({'Sc','Zr','Kl','Ca','Sp'})
legend({'[Ca^2^+]_c_y_t','[Ca^2^+]_e_x_t Steady State'},'location','bestoutside')
%legend({'[Ca^2^+]_c_y_t','[Ca^2^+]_e_x_t Steady','[Ca^2^+]_e_x_t Transient'},'location','bestoutside')
ylabel('MI(bit)')
%[hzcm pzcm]=ttest2(TTm{1},TTm{3})
%[hzcx pzcx]=ttest2(TTx{1},TTx{3})
%[hkcm pkcm]=ttest2(TTm{2},TTm{3})
%[hkcx pkcx]=ttest2(TTx{2},TTx{3})
[hicm picm]=ttest2([TTm{1:3}],[TTm{4:5}]);
[hicx picx]=ttest2([TTx{1:3}],[TTx{4:5}]);
disp([num2str(mean([TTm{1:3}])),'vs.',num2str(mean([TTm{4:5}])),',p=',num2str(picm)])
disp([num2str(mean([TTx{1:3}])),'vs.',num2str(mean([TTx{4:5}])),',p=',num2str(picx)])