%examples of single cell trajectories during trancient state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_4";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_4";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_4";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_4";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_4";
%TPLOT = cell(1,5);
cell_id = [35,221,56,105,57];
for i =[1:5]
    [TPLOT{i}]=SingleCellData(folder{i},1,[],cell_id(i));
end
%
%examples of single cell trajectories during steady state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_5";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_5";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_5";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_5";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_5";
SPLOT = cell(1,5);
cell_id = [244,331,294,155,175];
for i =[1:5]
    [SPLOT{i}]=SingleCellData(folder{i},1,[],cell_id(i));
end
%
%examples of single cell trajectories during no stress state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_1";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_1";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_1";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_1";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_26_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_26_1";
%NPLOT = cell(1,5);
cell_id = [144,131,194,68,108];
for i =[1:5]
    [NPLOT{i}]=SingleCellData(folder{i},1,[],cell_id(i));
end
%%
figure(2)
ind_t = 300;
for i = 1:5
    subplot(5,1,i)
    yyaxis left
    plot(TPLOT{i}{1}-ind_t,TPLOT{i}{2},'b-')%-- 2020-09-09 9:20 AM --%
    hold on
    plot(SPLOT{i}{1}+TPLOT{i}{1}(end)-ind_t,SPLOT{i}{2},'b-')
    plot(NPLOT{i}{1}-TPLOT{i}{1}(end)-ind_t,NPLOT{i}{2},'b-')

    ylim([-0.1 1])
    
    yyaxis right
    plot(TPLOT{i}{7}-ind_t,TPLOT{i}{8},'c.',...
        TPLOT{i}{7}-ind_t,TPLOT{i}{9},'k-','MarkerSize',0.5)
    hold on
    plot(SPLOT{i}{7}+TPLOT{i}{1}(end)-ind_t,SPLOT{i}{8},'c.',...
        SPLOT{i}{7}+TPLOT{i}{1}(end)-ind_t,SPLOT{i}{9},'k-','MarkerSize',0.5)
    plot(NPLOT{i}{7}-TPLOT{i}{1}(end)-ind_t,NPLOT{i}{8},'c.',...
        NPLOT{i}{7}-TPLOT{i}{1}(end)-ind_t,NPLOT{i}{9},'k-','MarkerSize',0.5)
    
    ylim([-0.1 0.5])
    
    ax = gca;%ax.YTick = [];
    if i~= 5, ax.XTickLabel = {''};end
    
    plot([0,0]-ind_t,[ax.YLim],'k--',[3600,3600]-ind_t,[ax.YLim],'k--')
    xlim([-3600 7200]-ind_t)
    if i == 3, yyaxis right; ylabel({'Nuclear Localization Score (AU)'});
               yyaxis left; ylabel('[Ca^2^+]_c_y_t (AU)');end
    if i == 5, xlabel('Time after calcium induction (sec)'); end
    yyaxis right; hold off
    yyaxis left; hold off
end
