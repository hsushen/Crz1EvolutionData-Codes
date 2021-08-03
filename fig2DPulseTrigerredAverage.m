%% ingroup

folders{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_5";
folders{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_5";
folders{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_13_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_13_3";
folders{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_17_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_17_3";
folders{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_5";
folders{6}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_4_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_04_3";
folders{7}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_7_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_07_3";
folders{8}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_5";
folders{9}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_20_2020\GCaMP3_CaIDR_dZF_GFP_2020_02_20_3";
folders{10}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_19_2020\GCaMP3_CaIDR_dZF_GFP_2020_03_19_3";
folders{11}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_5";
remove{1}=[9,11,15,30,45,59,64,67,69,70,74,80,92,94,103,111,118:120,138,142,149,151:153,161,165,169,182,193,205,207,208,218,223,225,248,251];
remove{2}=[31,45,97,106,121,155,182,266,296,243,354,389,437,444,449,];
remove{3}=[12,85,88,98,104,118,138,197,224,252,261,268,270];
remove{4}=[133];
remove{5}=[18,42,66,90,95,101,185,194,216,227,241,242,255,264,268,272,292,312,314,335,340,361,369,374,375,380];
remove{6}=[8,19,25,33,47,57,62,64,73,80,113,124,134,138,143,172,203,205,229,241,260,270];
remove{7}=[1,2,36,59,63,66,76,77,79,86,99,105,138,158,159,163,236,252,263,270,273];
remove{8}=[76,92,110,134,141,174,185,216];
remove{9}=[93];
remove{10}=[106,188];
remove{11}=[11,107,108,110,113,121,141,143,153,169,177,267];

for i = 1:5
    %pulse triggerrded averaging
    figure(1)
    num_cond=1;CICind=8;window=10*10;time = [-window:window]*6;
    if i == 1
        [Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n]=PulseTriggeredAverage(folders{1},num_cond,remove{1},CICind,window);
        subplot(1,5,i)
        Lcol = {[0 0 1],[1 0 0]};
        PlotPTA(time,Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n,Lcol)
        str = ['n=',num2str(sum(Ca_n(:,window+1)))];
        text(0,0.95,str,'Units','normalized')
        yyaxis left; ylabel('[Ca^2^+]_c_y_t (AU)')
    end
    if i == 2
        Ca_meanZr=nan(3,2*window+1);Ca_seZr=Ca_meanZr;Ca_nZr=Ca_meanZr;
        LS_meanZr=Ca_meanZr;LS_seZr=Ca_meanZr;LS_nZr=Ca_meanZr;
        for j = 2:4
            [Ca_meanZr(j-1,:),Ca_seZr(j-1,:),Ca_nZr(j-1,:),LS_meanZr(j-1,:),LS_seZr(j-1,:),LS_nZr(j-1,:)]=PulseTriggeredAverage(folders{j},num_cond,remove{j},CICind,window);
        end
        Ca_mean = sum(Ca_meanZr.*Ca_nZr)./sum(Ca_nZr);
        Ca_se = sum(Ca_seZr.*Ca_nZr.^1.5)./sum(Ca_nZr).^1.5;
        Ca_n = sum(Ca_nZr);
        LS_mean = sum(LS_meanZr.*LS_nZr)./sum(LS_nZr);
        LS_se = sum(LS_seZr.*LS_nZr.^1.5)./sum(LS_nZr).^1.5;
        LS_n = sum(LS_nZr);
        subplot(1,5,i)
        Lcol = {[0 0 1],[1 0 0]};
        PlotPTA(time,Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n,Lcol)
        str = ['n=',num2str(sum(Ca_n(:,window+1)))];
        text(0,0.95,str,'Units','normalized')
        
    end
     if i == 3
        Ca_meanKl=nan(3,2*window+1);Ca_seKl=Ca_meanKl;Ca_nKl=Ca_meanKl;
        LS_meanKl=Ca_meanKl;LS_seKl=Ca_meanKl;LS_nKl=Ca_meanKl;
        for j = 5:7
            [Ca_meanKl(j-4,:),Ca_seKl(j-4,:),Ca_nKl(j-4,:),LS_meanKl(j-4,:),LS_seKl(j-4,:),LS_nKl(j-4,:)]=PulseTriggeredAverage(folders{j},num_cond,remove{j},CICind,window);
        end
        Ca_mean = sum(Ca_meanKl.*Ca_nKl)./sum(Ca_nKl);
        Ca_se = sum(Ca_seKl.*Ca_nKl.^1.5)./sum(Ca_nKl).^1.5;
        Ca_n = sum(Ca_nKl);
        LS_mean = sum(LS_meanKl.*LS_nKl)./sum(LS_nKl);
        LS_se = sum(LS_seKl.*LS_nKl.^1.5)./sum(LS_nKl).^1.5;
        LS_n = sum(LS_nKl);
        subplot(1,5,i)
        Lcol = {[0 0 1],[1 0 0]};
        PlotPTA(time,Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n,Lcol)
        str = ['n=',num2str(sum(Ca_n(:,window+1)))];
        text(0,0.95,str,'Units','normalized')
        xlabel('Time after calcium bursts (sec)') 
     end
     if i == 4
        Ca_meanCa=nan(3,2*window+1);Ca_seCa=Ca_meanCa;Ca_nCa=Ca_meanCa;
        LS_meanCa=Ca_meanCa;LS_seCa=Ca_meanCa;LS_nCa=Ca_meanCa;
        for j = 8:10
            [Ca_meanCa(j-7,:),Ca_seCa(j-7,:),Ca_nCa(j-7,:),LS_meanCa(j-7,:),LS_seCa(j-7,:),LS_nCa(j-7,:)]=PulseTriggeredAverage(folders{j},num_cond,remove{j},CICind,window);
        end
        Ca_mean = sum(Ca_meanCa.*Ca_nCa)./sum(Ca_nCa);
        Ca_se = sum(Ca_seCa.*Ca_nCa.^1.5)./sum(Ca_nCa).^1.5;
        Ca_n = sum(Ca_nCa);
        LS_mean = sum(LS_meanCa.*LS_nCa)./sum(LS_nCa);
        LS_se = sum(LS_seCa.*LS_nCa.^1.5)./sum(LS_nCa).^1.5;
        LS_n = sum(LS_nCa);
        subplot(1,5,i)
        Lcol = {[0 0 1],[0 0 0]};
        PlotPTA(time,Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n,Lcol)
        str = ['n=',num2str(sum(Ca_n(:,window+1)))];
        text(0,0.95,str,'Units','normalized')
        
    end
    if i == 5
        [Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n]=PulseTriggeredAverage(folders{11},num_cond,remove{11},CICind,window);
        subplot(1,5,i)
        Lcol = {[0 0 1],[0 0 0]};
        PlotPTA(time,Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n,Lcol)
        str = ['n=',num2str(sum(Ca_n(:,window+1)))];
        text(0,0.95,str,'Units','normalized')
        yyaxis right; ylabel({'Nuclear Localization Score (AU)'})
    end
    
    
    
    
             
    
    
end
%%
function PlotPTA(time,Ca_mean,Ca_se,Ca_n,LS_mean,LS_se,LS_n,Lcol)
    
    yyaxis left
    Cax = [time,flip(time')'];
    Cay = ([Ca_mean + Ca_se*1.96, flip((Ca_mean- Ca_se*1.96)')']')';

    Cas = fill(Cax, Cay, Lcol{1});%,...
    Cas.EdgeColor ='none';
    Cas.FaceAlpha =0.15;
    hold on

    pCa=plot(time,Ca_mean,'-','Color',Lcol{1});
    hold off
    ylim([-0.1 0.75])
    yyaxis right
    LSx = [time,flip(time')'];
    LSy = ([LS_mean + LS_se*1.96, flip((LS_mean- LS_se*1.96)')']')';

    LSs = fill(LSx, LSy, Lcol{2});%,...
    LSs.EdgeColor ='none';
    LSs.FaceAlpha =0.15;
    hold on

    pLS=plot(time,LS_mean,'-','Color',Lcol{2});
    hold off
    ylim([0 0.25])
end


