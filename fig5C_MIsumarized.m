%Mutual information quantified between the two dynamics and between the stress and Crz1 dynamics 
%MI single cell trajectories during steady state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_5";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_5";
folder{3}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_13_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_13_3\';
folder{4}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_17_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_17_3\';
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_5";
folder{6}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_4_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_04_3\';
folder{7}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_7_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_07_3\';
folder{8}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_5";
folder{9}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_20_2020\GCaMP3_CaIDR_dZF_GFP_2020_02_20_3\';
folder{10}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_19_2020\GCaMP3_CaIDR_dZF_GFP_2020_03_19_3\';
folder{11}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_5";
CaMI = cell(1,11);
for i = 1:11
    cd(folder{i})
    CaMI{i} = load('Mi.mat');
end
% CaScCDS
folders{1}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_29_2020\GCaMP3_CaScCDSIDR_dZF_GFP_2020_02_29_3\';
folders{2}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_10_2020\GCaMP3_CaScCDSIDR_dZF_GFP_2020_03_10_3\';
folders{3}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_11_2020\GCaMP3_CaScCDSIDR_dZF_GFP_2020_03_11_3\';
CDSMI = cell(1,3);
for i = 1:3
    cd(folders{i})
    CDSMI{i}=load('MI.mat');
end
%

% Ca444Sc330
folders{1}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_27_2020\GCaMP3_Ca444Sc330IDR_dZF_GFP_2020_02_27_3\';
folders{2}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_12_2020\GCaMP3_Ca444Sc330IDR_dZF_GFP_2020_03_12_3\';
folders{3}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_13_2020\GCaMP3_Ca444Sc330IDR_dZF_GFP_2020_03_13_3\';
Ca444Sc330MI = cell(1,3);
for i = 1:3
    cd(folders{i})
    Ca444Sc330MI{i}=load('MI.mat');
end
%

% Ca451Sc337
folders{1}='C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_28_2020\GCaMP3_Ca451Sc337IDR_dZF_GFP_2020_02_28_3\';
Ca451Sc337MI = cell(1,3);
for i = 1
    cd(folders{i})
    Ca451Sc337MI{i}=load('MI.mat');
end

%
m = nan(2,1); errm = m;
%% load MI of GP value
cd('C:\Users\Ian\Documents\MATLAB\Crz1\MI quantification\SteadyState')
SMI = load('MIReporter.mat');
SMIch = load('MIReporterChimera.mat');
x = [
    mean([SMI.ScMI,SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3,SMI.KlMI1,SMI.KlMI2,SMI.KlMI3]);
    mean([SMI.CaMI1,SMI.CaMI2,SMI.CaMI3,SMI.SpMI]);
    ];
errx = [ std([SMI.ScMI,SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3,SMI.KlMI1,SMI.KlMI2,SMI.KlMI3]);
    std([SMI.CaMI1,SMI.CaMI2,SMI.CaMI3,SMI.SpMI])]/3^0.5;
TTx = { [SMI.ScMI,SMI.ZrMI1,SMI.ZrMI2,SMI.ZrMI3,SMI.KlMI1,SMI.KlMI2,SMI.KlMI3];
    [SMI.CaMI1,SMI.CaMI2,SMI.CaMI3,SMI.SpMI]};
xch = [
    mean([SMIch.CaSc4MI1,SMIch.CaSc4MI2,SMIch.CaSc4MI3]);
    mean([SMIch.CaSc2MI1,SMIch.CaSc2MI2,SMIch.CaSc2MI3]);
    SMIch.CaSc3MI;
        ];
errxch = [
        std([SMIch.CaSc4MI1,SMIch.CaSc4MI2,SMIch.CaSc4MI3]);
        std([SMIch.CaSc2MI1,SMIch.CaSc2MI2,SMIch.CaSc2MI3]);
        0;
        ]/3^0.5;
TTxch = {
        [SMIch.CaSc4MI1,SMIch.CaSc4MI2,SMIch.CaSc4MI3];
        [SMIch.CaSc2MI1,SMIch.CaSc2MI2,SMIch.CaSc2MI3];
        [SMIch.CaSc3MI];};
xall = [x;xch];
errxall = [errx;errxch];
%%
TTm = cell(1,2); 
ami = zeros(11,30);
for i = 1:11
    ami(i,:) = mean(CaMI{i}.I);
end
%
    m(1) = mean(ami(1:7,1));
    m(2) = mean(ami(8:11,1));
    errm(1)=std(ami(1:7,1),1);
    errm(2)=std(ami(8:11,1),1);
    TTm{1}=[ami(1:7,1)];
    TTm{2}=[ami(8:11,1)];
%
mch = zeros(3,1); errmch = mch;
TTmch = cell(1,3);
ami2 = zeros(3,2);
for i = 1:3
    amicds = mean(CDSMI{i}.I);
    ami444 = mean(Ca444Sc330MI{i}.I);
    ami2(i,:) = [amicds(1) ami444(1)];
end
mch(1:2)=mean(ami2);errmch(1:2)=std(ami2,1);
ami451 = mean(Ca451Sc337MI{1}.I); 
mch(3) = ami451(1);

TTmch{1}= ami2(:,1);TTmch{2}= ami2(:,2);TTmch{3} = mch(3);
%
mall = [m;mch];
figure(4)
bar([mall,xall])
hold on
dif = 0.15;
errmall = [errm; errmch];n = [7,4,3,3,1]'.^0.5;
er1 = errorbar([1:length(mall)]'-dif,mall,errmall./n*1.96);
er1.Color = [0 0 0];
er1.LineStyle = 'none';
er1 = errorbar([1:length(xall)]'+dif,xall,errxall./n*1.96);
er1.Color = [0 0 0];
er1.LineStyle = 'none';
hold off

xticklabels({'Saccharomyces','Outgroup','Ca^H^i^g^h','Ca^P^x^I^x^I^T^+^N^L^S','Ca^N^L^S'})
ylabel('MI(bit)')
% m
[hOP pOP]=ttest2([TTm{2}],[TTmch{1}]);
[hIP pIP]=ttest2([TTm{1}],[TTmch{1}]);

[hON pON]=ttest2([TTm{2}],[TTmch{2}]);
[hIN pIN]=ttest2([TTm{1}],[TTmch{2}]);
hold on 
%plot([1,3]'-dif,[1,1]*0.6,'k-')
%plot([1,4]'-dif,[1,1]*0.65,'k-')
plot([2,3]'-dif,[1,1]*-0.1,'k-')
plot([2,4]'-dif,[1,1]*-0.15,'k-')
hold off

%text(2-dif,0.6,['N.S'],'VerticalAlignment','bottom','HorizontalAlignment','center')
%text(2.5-dif,0.65,['N.S'],'VerticalAlignment','bottom','HorizontalAlignment','center')
text(2.5-dif,-0.1,['p=',num2str(round(pOP*10^5)/10^5)],'VerticalAlignment','bottom','HorizontalAlignment','center')
text(3-dif,-0.15,['p=',num2str(round(pON*10^5)/10^5)],'VerticalAlignment','bottom','HorizontalAlignment','center')

[hIO pIO]=ttest2([TTm{2}],[TTm{1}])

% x
[hOP pOP]=ttest2([TTx{2}],[TTxch{1}]);
[hIP pIP]=ttest2([TTx{1}],[TTxch{1}]);

[hON pON]=ttest2([TTx{2}],[TTxch{2}]);
[hIN pIN]=ttest2([TTx{1}],[TTxch{2}]);
hold on 
%plot([1,3]'+dif,[1,1]*1,'k-')
%plot([1,4]'+dif,[1,1]*1.05,'k-')
%plot([2,3]'+dif,[1,1]*0.8,'k-')
%plot([2,4]'+dif,[1,1]*0.85,'k-')
hold off

% text(2+dif,1,['p=',num2str(round(pIP*10^5)/10^5)],'VerticalAlignment','bottom','HorizontalAlignment','center')
% text(2.5+dif,1.05,['p=',num2str(round(pIN*10^5)/10^5)],'VerticalAlignment','bottom','HorizontalAlignment','center')
% text(2.5+dif,0.8,['p=',num2str(round(pOP*10^5)/10^5)],'VerticalAlignment','bottom','HorizontalAlignment','center')
% text(3+dif,0.85,['p=',num2str(round(pON*10^5)/10^5)],'VerticalAlignment','bottom','HorizontalAlignment','center')

%%
function [fI,fIn]= ThreeR(MI,col,eb,name)
%subplot(1,2,2)
    x = [1:30]*0.5;
    y = [mean(MI{1}.I);mean(MI{2}.I);mean(MI{3}.I)];
    
    xfit = y; for i = 1:size(y,2), xfit(:,i)= x(i);end
    fI = fit(xfit(:),y(:),'exp1');
   % plot(xfit(:),y(:),[col{1},'.'],'MarkerSize',1)
    %ylim([-0.1 1])    
    if eb==1
        err = std([mean(MI{1}.I);mean(MI{2}.I);mean(MI{3}.I)],1)/3^0.5*1.96;
        hold on
       %errorbar(x,mean(y),err,[col,'-'])
        hold off
    end
    
 %   subplot(1,2,1)
 hold on
    xn = [1:30]*-0.5;
    yn = [mean(MI{1}.In);mean(MI{2}.In);mean(MI{3}.In)];
    xfitn = yn; for i = 1:size(yn,2), xfitn(:,i)= xn(i);end
    fIn = fit(xfitn(:),yn(:),'exp1');
    %plot(fIn,xfit(:),y(:),[col,'.'])
   % plot(xfitn(:),yn(:),[col{1},'.'],'MarkerSize',1)
    %ylim([-0.1 0.8])
    if eb==1
        errn = std([mean(MI{1}.In);mean(MI{2}.In);mean(MI{3}.In)],1);%/3^0.5*1.96;
        hold on
        %errorbar([flip(xn),x],mean([flip(yn),y]),[flip(errn),err],[col,'-'])
        errorbar([flip(xn),x],[flip(mean(yn)),mean(y)],[flip(errn),err],[col{1},col{2}],'DisplayName',name)
        hold off
    end
    
    hold on
   % plot([flip(xn),x],[flip(mean(yn)),mean(y)],[col{1},col{2}])
    %plot([-15,15],[0,0],'k--')
    %xlim([-15,15])
    %xlim([-5,5])
    hold off
    
end
%%
function [fI,fIn]= OneR(MI,col,eb,name)
%subplot(1,2,2)
    x = [1:30]*0.5;
    y = [mean(MI{1}.I)];
    
    xfit = y; for i = 1:size(y,2), xfit(:,i)= x(i);end
    fI = fit(xfit(:),y(:),'exp1');
    %plot(xfit(:),y(:),[col{1},'.'],'MarkerSize',1)
    %ylim([-0.1 1])    
    if eb==1
        %err = std([mean(MI{1}.I);mean(MI{2}.I);mean(MI{3}.I)],1)/3^0.5*1.96;
        hold on
       %errorbar(x,mean(y),err,[col,'-'])
        hold off
    end
    
 %   subplot(1,2,1)
 hold on
    xn = [1:30]*-0.5;
    yn = [mean(MI{1}.In)];
    xfitn = yn; for i = 1:size(yn,2), xfitn(:,i)= xn(i);end
    fIn = fit(xfitn(:),yn(:),'exp1');
    %plot(fIn,xfit(:),y(:),[col,'.'])
    %plot(xfitn(:),yn(:),[col{1},'.'],'MarkerSize',1)
    %ylim([-0.1 0.8])
    if eb==1
        %errn = std([mean(MI{1}.In);mean(MI{2}.In);mean(MI{3}.In)],1)/3^0.5*1.96;
        hold on
        %errorbar([flip(xn),x],mean([flip(yn),y]),[flip(errn),err],[col,'-'])
        %errorbar([flip(xn),x],[flip(mean(yn)),mean(y)],[flip(errn),err],[col,'-'])
        hold off
    end
    
    hold on
    plot([flip(xn),x],[flip(yn),y],[col{1},col{2}],'DisplayName',name)
    %plot([-15,15],[0,0],'k--')
    %xlim([-15,15])
    %xlim([-5,5])
    hold off
    
end