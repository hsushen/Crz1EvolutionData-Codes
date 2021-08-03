% populational trajectories of every orthologous reporter no stress
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_1";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_1";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_13_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_13_1";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_17_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_17_1";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_1";
folder{6}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_4_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_04_1";
folder{7}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_7_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_07_1";
folder{8}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_1";
folder{9}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_20_2020\GCaMP3_CaIDR_dZF_GFP_2020_02_20_1";
folder{10}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_19_2020\GCaMP3_CaIDR_dZF_GFP_2020_03_19_1";
folder{11}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_1";
[Nt1,Nmean_l,Nse_l,Nn]=ALLExp(folder);

% populational trajectories of every orthologous reporter transient state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_4";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_4";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_13_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_13_2";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_17_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_17_2";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_4";
folder{6}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_4_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_04_2";
folder{7}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_7_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_07_2";
folder{8}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_4";
folder{9}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_20_2020\GCaMP3_CaIDR_dZF_GFP_2020_02_20_2";
folder{10}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_19_2020\GCaMP3_CaIDR_dZF_GFP_2020_03_19_2";
folder{11}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_4";
[Tt1,Tmean_l,Tse_l,Tn]=ALLExp(folder);

% populational trajectories of every orthologous reporter steady state
addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')
folder{1}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_17_2019\GCaMP3_ScIDR_dZF_GFP_2019_05_17_5";
folder{2}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_7_2019\GCaMP3_ZrIDR_dZF_GFP_2019_05_07_5";
folder{3}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_13_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_13_3";
folder{4}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_17_2020\GCaMP3_ZrIDR_dZF_GFP_2020_02_17_3";
folder{5}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_5_10_2019\GCaMP3_KlIDR_dZF_GFP_2019_05_10_5";
folder{6}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_4_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_04_3";
folder{7}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_7_2020\GCaMP3_KlIDR_dZF_GFP_2020_02_07_3";
folder{8}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_8_26_2019\GCaMP3_CaIDR_dZF_GFP_2019_08_26_5";
folder{9}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_20_2020\GCaMP3_CaIDR_dZF_GFP_2020_02_20_3";
folder{10}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_19_2020\GCaMP3_CaIDR_dZF_GFP_2020_03_19_3";
folder{11}="C:\Users\Ian\Documents\MATLAB\Crz1\crz1_9_3_2019\GCaMP3_SpIDR_dZF_GFP_2019_09_03_5";
[St1,Smean_l,Sse_l,Sn]=ALLExp(folder);



%
figure(2)
Lcol = {[1 0 0],[1 0 1],[0 1 0],[0 1 1],[0 0 1]};
leg = {'Sc','Zr','Kl','Ca','Sp'};
in_t = 300;
%y_l = [mean_l + se_l, flip((mean_l- se_l)')'];
%Transient state
%subplot(1,2,1)
Tx1 = [Tt1,flip(Tt1')']-in_t;
Tyl = ([Tmean_l + Tse_l, flip((Tmean_l- Tse_l)')']')';
for i = 1:5
    s_l = fill(Tx1(i,:), Tyl(i,:), Lcol{i});%,...
    s_l.EdgeColor ='none';
    s_l.FaceAlpha =0.15;
hold on
end

for i = 1:5
p{i} = plot(Tt1(i,:)'-in_t,Tmean_l(i,:)','Color',Lcol{i});
%plot(t1(i,:)',mean_l(i,:)','Color',Lcol{i})

end
ax = gca;
%xlim([0 3600])

%Steady state

Sx1 = [St1+Tt1(i,end),flip((St1+Tt1(i,end))')']-in_t;
Syl = ([Smean_l + Sse_l*1.96, flip((Smean_l- Sse_l*1.96)')']')';
for i = 1:5
    s_l = fill(Sx1(i,:), Syl(i,:), Lcol{i});%,...
    s_l.EdgeColor ='none';
    s_l.FaceAlpha =0.15;
hold on
end

for i = 1:5
p{i} = plot(St1(i,:)'+Tt1(i,end)-in_t,Smean_l(i,:)','Color',Lcol{i});
%plot(t1(i,:)',mean_l(i,:)','Color',Lcol{i})

end
%No Stress state

Nx1 = [Nt1-Tt1(i,end),flip((Nt1-Tt1(i,end))')']-in_t;
Nyl = ([Nmean_l + Nse_l, flip((Nmean_l- Nse_l)')']')';
for i = 1:5
    s_l = fill(Nx1(i,:), Nyl(i,:), Lcol{i});%,...
    s_l.EdgeColor ='none';
    s_l.FaceAlpha =0.15;
hold on
end

for i = 1:5
p{i} = plot(Nt1(i,:)'-Tt1(i,end)-in_t,Nmean_l(i,:)','Color',Lcol{i});
%plot(t1(i,:)',mean_l(i,:)','Color',Lcol{i})

end

plot([0,0]-in_t,[ax.YLim],'k--',[3600,3600]-in_t,[ax.YLim],'k--')

hold off
legend([p{:}],leg{:},'location','bestoutside');
ylim(ax.YLim)
xlabel('Time after calcium induction (sec)')
ylabel({'Localization';'Score (AU)'})
xlim([-3600 7200]-in_t)
%%
function [t1,mean_l,se_l,n]=ALLExp(folder)
    t1 = nan(5,601);mean_l = t1; se_l = t1; n = nan(5,1);
    for i = 1:5
        if i == 1
        [~,~,~,t1(i,:),mean_l(i,:),se_l(i,:),n(i)]=TransientData(folder{1},1,[]);
        end
        if i == 2
            t1Zr = nan(3,601);mean_lZr = t1Zr; se_lZr = t1Zr; nZr = nan(3,1);
            for j = 2:4
                [~,~,~,t1Zr(j-1,:),mean_lZr(j-1,:),se_lZr(j-1,:),nZr(j-1)]=TransientData(folder{j},1,[]);
            end
            t1(i,:) = sum(t1Zr.*nZr)/sum(nZr); 
            mean_l(i,:) = sum(mean_lZr.*nZr)/sum(nZr);
            se_l(i,:) = sum(se_lZr.*nZr.^1.5)/sum(nZr)^1.5;
        end
        if i == 3
            t1Kl = nan(3,601);mean_lKl = t1Kl; se_lKl = t1Kl; nKl = nan(3,1);
            for j = 5:7
                [~,~,~,t1Kl(j-4,:),mean_lKl(j-4,:),se_lKl(j-4,:),nKl(j-4)]=TransientData(folder{j},1,[]);
            end
            t1(i,:) = sum(t1Kl.*nKl)/sum(nKl); 
            mean_l(i,:) = sum(mean_lKl.*nKl)/sum(nKl);
            se_l(i,:) = sum(se_lKl.*nKl.^1.5)/sum(nKl)^1.5;
        end
        if i == 4
            t1Ca = nan(3,601);mean_lCa = t1Ca; se_lCa = t1Ca; nCa = nan(3,1);
            for j = 8:10
                [~,~,~,t1Ca(j-7,:),mean_lCa(j-7,:),se_lCa(j-7,:),nCa(j-7)]=TransientData(folder{j},1,[]);
            end
            t1(i,:) = sum(t1Ca.*nCa)/sum(nCa); 
            mean_l(i,:) = sum(mean_lCa.*nCa)/sum(nCa);
            se_l(i,:) = sum(se_lCa.*nCa.^1.5)/sum(nCa)^1.5;
        end
        if i == 5
            [~,~,~,t1(i,:),mean_l(i,:),se_l(i,:),n(i)]=TransientData(folder{11},1,[]);
        end
    end
end