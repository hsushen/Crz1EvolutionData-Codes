%examples of single cell trajectories during no stress state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_27_2020\GCaMP3_Ca444Sc330IDR_dZF_GFP_2020_02_27_1";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_28_2020\GCaMP3_Ca451Sc337IDR_dZF_GFP_2020_02_28_1";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_29_2020\GCaMP3_CaScCDSIDR_dZF_GFP_2020_02_29_1";
NPLOT = cell(1,3);
cell_id = [103,120,160];
for i =[1:3]
    [NPLOT{i}]=SingleCellData(folder{i},1,[],cell_id(i));
end
%
%examples of single cell trajectories during trancient state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_27_2020\GCaMP3_Ca444Sc330IDR_dZF_GFP_2020_02_27_2";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_28_2020\GCaMP3_Ca451Sc337IDR_dZF_GFP_2020_02_28_2";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_29_2020\GCaMP3_CaScCDSIDR_dZF_GFP_2020_02_29_2";
TPLOT = cell(1,3);
cell_id = [82,121,151];
for i =[1:3]
    [TPLOT{i}]=SingleCellData(folder{i},1,[],cell_id(i));
end
%
%examples of single cell trajectories during steady state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_27_2020\GCaMP3_Ca444Sc330IDR_dZF_GFP_2020_02_27_3";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_28_2020\GCaMP3_Ca451Sc337IDR_dZF_GFP_2020_02_28_3";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_29_2020\GCaMP3_CaScCDSIDR_dZF_GFP_2020_02_29_3";
SPLOT = cell(1,3);
cell_id = [195,196,132];
for i =[1:3]
    [SPLOT{i}]=SingleCellData(folder{i},1,[],cell_id(i));
end
%

%%
figure(2)
ind_t = 300;
for i = 1:3
    subplot(3,1,i)
    yyaxis left
    plot(TPLOT{i}{1}-ind_t,TPLOT{i}{2},'b-')
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
    if i~= 3, ax.XTickLabel = {''};end
    
    plot([0,0]-ind_t,[ax.YLim],'k--',[3600,3600]-ind_t,[ax.YLim],'k--')
    xlim([-3600 7200]-ind_t)
    if i == 2, yyaxis right; ylabel({'Nuclear Localization Score (AU)'});
               yyaxis left; ylabel('[Ca^2^+]_c_y_t (AU)');end
    if i == 3, xlabel('Time after calcium induction (sec)'); end
    yyaxis right; hold off
    yyaxis left; hold off
end
