% normalized s for fitness strains
cd('C:\Users\Ian\Google Drive\Moses Lab\pulsing\CRZ1')
T = readtable('Crz1 analysis.xlsx','Sheet','Fitness');
Ts = convertvars(T(:,[1,3,4,6:11]),{'IDR','Stress','RBTag'},'categorical');
Day = {'2020-03-20'};
Strain = {'CDS','Ca'};
Strain = flip(Strain);
Cond = {'0.2M CaCl2','No stress'};
s = nan(length(Strain),2);se = s;
k = 1;
for i = 1:length(Strain)
    for j = 1:2
        pick = Ts.Day==Day&Ts.IDR==Strain{i}&Ts.Stress == Cond{j};%&Ts.RBTag ~='Red';
        %s(i,j) = mean(Ts.NormalizedS(pick))*-1;
        s(i,j) = mean(Ts.RBs(pick));
        N(k) = sum(Ts.n(pick));
        se(i,j) = (sum(Ts.SquareSumDeviation(pick))/(N(k)-1))^0.5/N(k)^0.5*1.96;%1.96 standard error
        k=k+1;
    end
end
%
figure(6)
bar(s)
hold on 
dif = 0.15;
er1 = errorbar([1:length(Strain)]'-dif,s(:,1),se(:,1));%,'horizontal');
er1.Color = [0 0 0];
er1.LineStyle = 'none';
er2 = errorbar([1:length(Strain)]'+dif,s(:,2),se(:,2));%,'horizontal');
er2.Color = [0 0 0];
er2.LineStyle = 'none';
hold off
xticklabels({'Ca','Ca^C^D^S'})
ylabel('s')

title(['n =', num2str(N)])

filename='C:\Users\Ian\Documents\MATLAB\Crz1\fitness_assay\setup 15-2\s.mat';
load(filename,'CDSO1sNS','CaOsNS','CDSO1sWS','CaOsWS')

[h1,p1]=ttest2(CDSO1sWS,CaOsWS)
[h2,p2]=ttest2(CDSO1sNS,CaOsNS);
hold on 
plot([1:length(Strain)]'-dif,ones(size(Strain))*0.07,'k-')
plot([1:length(Strain)]'+dif,ones(size(Strain))*-0.11,'k-')
hold off
text(1.5-dif,0.07,['p<',num2str(0.005)],'VerticalAlignment','bottom','HorizontalAlignment','center')
text(1.5+dif,-0.11,['N.S'],'VerticalAlignment','bottom','HorizontalAlignment','center')


legend(Cond,'Location','best')