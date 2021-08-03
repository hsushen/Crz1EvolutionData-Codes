% normalized s for fitness strains
cd('C:\Users\Ian\Google Drive\Moses Lab\pulsing\CRZ1')
T = readtable('Crz1 analysis.xlsx','Sheet','Fitness');
Ts = convertvars(T(:,[1,3,4,6:11]),{'IDR','Stress','RBTag'},'categorical');
%Strain = {'Sc','Zr','Kl','Ca','Sp','KO','mSRR','3A'};
Strain = {'Sc','Zr','Kl','Ca','Sp','mSRR','KO'};
Strain = Strain;
Cond = {'0.2M CaCl2','No stress'};
Day = {'2020-03-20'};
s = nan(length(Strain),2);se = s;
k = 1;
for i = 1:length(Strain)
    for j = 1:2
        pick = Ts.IDR==Strain{i}&Ts.Stress == Cond{j}&Ts.Day~=Day;%&Ts.RBTag ~='Red';
        %s(i,j) = mean(Ts.NormalizedS(pick))*-1;
        s(i,j) = mean(Ts.RBs(pick));
        N(k) = sum(Ts.n(pick));
        se(i,j) = (sum(Ts.SquareSumDeviation(pick))/(N(k)-1))^0.5/N(k)^0.5*1.96;
        k = k+1;
    end
end
%
figure(5)
bar(s)
hold on 
dif = 0.15;
er1 = errorbar([1:length(Strain)]'-dif,s(:,1),se(:,1));
er1.Color = [0 0 0];
er1.LineStyle = 'none';
er2 = errorbar([1:length(Strain)]'+dif,s(:,2),se(:,2));
er2.Color = [0 0 0];
er2.LineStyle = 'none';
hold off
xticklabels(Strain)
ylabel('s')

title(['n =', num2str(N)])
Tzc = (s(2,1)-s(4,1))/((se(2,1)/1.96*N(3)^0.5)^2/N(3)+(se(4,1)/1.96*N(7)^0.5)^2/N(7))^0.5;
Tkc = (s(3,1)-s(4,1))/((se(3,1)/1.96*N(5)^0.5)^2/N(5)+(se(4,1)/1.96*N(7)^0.5)^2/N(7))^0.5;
pzc = 1-tcdf(Tzc,N(3)+N(7)-2);
pkc = 1-tcdf(Tkc,N(5)+N(7)-2);
hold on 
l1 = -0.035;l2 = -0.045;
plot([2 4]'-dif,[1 1]*l1,'k-')
plot([3 4]'-dif,[1 1]*l2,'k-')
hold off
text(3-dif,l1,['p<',num2str(ceil(pzc*10^6)/10^6)],'VerticalAlignment','bottom','HorizontalAlignment','center')
text(3.5-dif,l2,['p<',num2str(ceil(pkc*10^4)/10^4)],'VerticalAlignment','bottom','HorizontalAlignment','center')
legend(Cond,'Location','best')