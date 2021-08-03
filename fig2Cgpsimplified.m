%Quantitative features of single cell trajectories
%%
%GP values of single cell trajectories during steady state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_5";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_5";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_5";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_5";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_5";
SGP = cell(1,5);
for i = 1:5
    cd(folder{i})
    SGP{i} = load('GP.mat');
end
%% violin plot
figure(1)
Lcol = [1 0 0;1 0 0;1 0 0;;0 0 0;0 0 0];
v = [1,2];logON = [1,1];xl = [3 8.2];yl = [-6 -1];
DynamicSpaceViolin(SGP,Lcol,v,logON,xl,yl)

%%
function DynamicSpaceViolin(SGP,Lcol,v,logON,xl,yl)
    n = nan(1,5);xdata = cell(1,5);ydata = xdata
    for i = 1:5
        xdata{i} = SGP{i}.V(v(1),:);if logON(1)==1, xdata{i} = log(xdata{i});end
        ydata{i} = SGP{i}.V(v(2),:);if logON(2)==1, ydata{i} = log(ydata{i});end
        xdata{i}(ydata{i}<-6)=[];ydata{i}(ydata{i}<-6)=[];n(i)=length(ydata{i});
    end
    Ya =nan(max(n),5); Yl = Ya;
    for i = 1:5
        Yl(1:length(xdata{i}),i)=xdata{i};
        Ya(1:length(ydata{i}),i)=ydata{i};
    end
    
    
    %xpoint = (1:5); 
    
    subplot(2,1,1)
    [h,L,MX,MED]=violin(Yl,'facecolor',Lcol);%,'x',xpoint);
    xticklabels({'Sc','Zr','Kl','Ca','Sp'})
    ylabel('log(l)')
    subplot(2,1,2)
    [h,L,MX,MED]=violin(Ya,'facecolor',Lcol);%,'x',xpoint);
    xticklabels({'Sc','Zr','Kl','Ca','Sp'})
    ylabel('log(a)')
%     leg = {['Sc',',n=',num2str(n(1))],...
%         ['Zr',',n=',num2str(n(2))],...
%         ['Kl',',n=',num2str(n(3))],...
%         ['Ca',',n=',num2str(n(4))],...
%         ['Sp',',n=',num2str(n(5))]};
%     legend([e{:}],leg{:},'location','northoutside');
%     
%     
%     xlim(xl)
%     ylim(yl)
end