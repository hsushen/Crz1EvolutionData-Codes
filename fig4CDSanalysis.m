folder = 'C:/Users/Ian/Documents/MATLAB/Crz1/orthologues';
cd(folder)
set(0,'DefaultFigureWindowStyle','normal') 

addpath('C:\Users\Ian\Documents\MATLAB\Crz1\Cell_trace_and_quantification')

%% CDS affinity change 3 (confident alignment;Sach and Ca clade;yeast-only PSSM)
clade = [78,55];%78:Sach clade; 55: Ca clade; 
motif_length = [1:12];
figure(2)
EdStest = 0;selecON=0;W=0;
%disp(['No Selection'])
[S3,ance_all3,DIS3,VdF,Z3,ID3,EdS3,nullEdS3]=CDSChange3(folder,clade,motif_length,EdStest,selecON,W);

if EdStest==1&&selecON ==0, save(['W_v3ns.mat'],'EdS3','nullEdS3'); end
Species = ID3(1:40,1);
if 1 == 0 %simulation for selection
EdStest = 1;selecON=1;
for W=-5.608%-6:0
    disp(['W=',num2str(W)])
    [~,~,~,~,~,~,EdS3,nullEdS3]=CDSChange3(folder,clade,motif_length,EdStest,selecON,W);

    if selecON ==1, save(['W_v3',num2str(W),'.mat'],'EdS3','nullEdS3'); end
end
end
Species = ID3(1:40,1);
% No Selection
% z test on E[dS] of each clade
% E[dS]Sac=0.58003,p =5.8354e-06,z =4.5323
% E[dS]Can=0.17516,p =0.0024727,z =3.0267
% E[dS]null=-1.3089+-5.5405 for Sacch,-1.2529+-5.5485 for Ca
% W=-6
% z test on E[dS] of each clade
% E[dS]Sac=0.58003,p =0.0063752,z =2.7278
% E[dS]Can=0.17516,p =0.19611,z =1.2927
% E[dS]null=-0.48654+-5.1977 for Sacch,-0.39295+-5.1681 for Ca
% W=0
% z test on E[dS] of each clade
% E[dS]Sac=0.58003,p =0.04867,z =1.9715
% E[dS]Can=0.17516,p =0.47133,z =0.72032
% E[dS]null=-0.033197+-4.1349 for Sacch,-0.081943+-4.1975 for Ca


%% plot CDS across time
Sc = 17;Zr = 29; Kl = 32; Ca = 2;  rest = setdiff(1:40,[Sc,Zr,Kl,Ca]);
Strains = [Sc,Zr,Kl,Ca,rest];
StX = cell(size(Strains));StY = StX;AAC = StX;
X = [DIS3(end)];Y = [ID3{end,4}];
for st = 1:length(Strains)
    ance = ance_all3{Strains(st)};
    X = [DIS3(end)];Y = [ID3{end,4}];
    aac = [];ch = [];
    for i = 2:length(ance)

        if isempty(ID3{ance(i),2})
            Xs = [X(end), DIS3(ance(i))+X(end)];
            Ys = [Y(end), Y(end)];
            X = [X, Xs];
            Y = [Y, Ys];
            aac = [aac, {''}];
        else 
            steps = DIS3(ance(i))/(length(ID3{ance(i),2})+1);
            changes = ID3{ance(i),2};
            aac_all = ID3{ance(i),3};
            for j = 1:length(changes)
                Xs = [X(end), steps + X(end)];
                Ys = [changes(j)+Y(end),changes(j)+Y(end)];
                X = [X, Xs];
                Y = [Y, Ys];
                
            end
            aac = [aac, aac_all];
            ch = [ch,changes];
        end
        
    end
    StX{st}=X;
    StY{st}=Y;
    AAC{st}= aac;
    disp(['E[ch]=',num2str(mean(ch))])
end
%
figure(3)
plot(StX{1},StY{1},'r-')
hold on
plot(StX{2},StY{2},'r--')
plot(StX{3},StY{3},'r--')
plot(StX{4},StY{4},'k--')

for i = [1:15]+4
    plot(StX{i},StY{i},'k:')
end
for i = [1:21]+19
    plot(StX{i},StY{i},'r:')
end

% label aa change
for st = 1:4;%[1:21]+19
    ance = ance_all3{Strains(st)};
    txX = StX{st}; txX = reshape(txX(1:end-1),2,[]);txX = txX(1,:);
    txY = StY{st}; txY = reshape(txY(2:end),2,[]);txY = txY(1,:);
    text(txX,txY+1,AAC{st})
end
hold off
xlabel('Time(AU)')
ylabel({'Predicted calcineurin'; 'docking strength (AU)'})
%ylim([-4 12])
%%
figure(4)
bin = 2.25;
histogram([ID3{17:40,4}],-4:bin:12,'facecolor','r','orientation','horizontal','normalization','probability')
hold on 
histogram([ID3{1:16,4}],-4:bin:12,'facecolor','k','orientation','horizontal','normalization','probability')
hold off
ylim([-4 12])
ylabel({'Predicted calcineurin'; 'docking strength (AU)'})
xlabel('Probability')
%% Z-test of E[dS] on pulsing IDRs
ance_pulse = unique([ance_all3{Sc},ance_all3{Zr},ance_all3{Kl}]);
dS_pulse = [ID3{ance_pulse,2}];
[h,p,~,z]=ztest(dS_pulse,0,VdF(4)^0.5);
disp(['E[dS]pulse=',num2str(mean(dS_pulse)),',p =',num2str(p),',z =',num2str(z)])

%% CDS affinity correlation between prediction and experimental data
figure(5),
[mdl] = CDScorr([1:12],0)
ylabel({'Predicted calcineurin'; 'docking strength (AU)'})
xlabel({'Measured calcineurin'; 'docking strength (AU)'})
%ylim([-4 12])
%% Fisher's exact test
dSnull = load('W_v3ns.mat','nullEdS3');
dS = load('W_v3ns.mat','EdS3');
edg = [-8,-1,1,8];
bin0 = 2;bin1 = 3;
figure(4),
ST = {'Sc','Ca'};
for st = 1:2
    subplot(2,1,st)
    
histogram(dSnull.nullEdS3{st},edg,'normalization','probability')
hold on
histogram(dS.EdS3{st},edg,'normalization','probability')
hold off
title([ST{st},' clade binned dS'])
xlabel('dS')
ylabel('Probability')
bindS = histcounts(dS.EdS3{st},edg);%,'normalization','probability')
bindSnu = histcounts(dSnull.nullEdS3{st},edg);%,'normalization','probability')
x = table([bindSnu(bin0);bindS(bin0) ],[bindSnu(bin1);bindS(bin1) ],'VariableNames',{'dS0','dS8'},'RowNames',{'NoSelection','Sacch'});
[h,p] = fishertest(x);
disp([ST{st},'|p=',num2str(p)])
end
%Sc|p=0.023817
%Ca|p=0.67835
%% CDS affinity 3 (Sc and Ca only, 40 species)
fixregion = 1;
% define region relative to ScCrz1
%scregion = [331:430];%NES:186-279;CDS:331-2;CaSc:392;ZF:568
%scregion = [330:429];
binsize = 4;
psite = 0;
motif_length = [1:12];


scregion = cell(1,2);
scregion{1} = [1:568];
%scregion{2} = [1:(230-53)];
%scregion{3} = [231:330]-53;
scregion{2} = [331:430]-53;
%scregion{5} = [431:530]-53;
HMall = nan(2,40);
for i = 1:2
    figure(6),
[CDSA,ID]=cdsaffinity3(folder,fixregion,scregion{i},binsize,psite,motif_length);
HMall(i,:) = CDSA;
Species = flip({ID{1:40}})';
T = table(Species,CDSA);
writetable(T,['Crz1CDSA',num2str(scregion{i}(1)),'to',num2str(scregion{i}(end)),'_v4.csv'],'Delimiter',',','QuoteStrings',true)
end
%% combine affinity from aligned sequence
HMall(3,:) = flip([ID3{1:size(HMall,2),4}]);
%% heatmap of every region
%hmscale = [0 floor(max(HMall(:)))];
hmscale = [0 7];
figure(7)
HMallsc = HMall; 
HMallsc(HMallsc<hmscale(1)) = hmscale(1); 
HMallsc(HMallsc>hmscale(2)) = hmscale(2);
b=bar3(flip(HMallsc-min(HMallsc(:)),2));%ones(size(Headers)),flip(1:length(Headers)),
for i = 1:40
zdata = b(i).ZData;
b(i).CData = zdata;%(zdata-min(zdata))./(max(zdata)-min(zdata))+1;
b(i).FaceColor = 'interp';b(i).EdgeColor = 'none';
end
view(2)

ax = gca; ax.XColor = 'none';
ax.YColor = 'none';ax.Color = 'none';ax.XGrid = 'off'; ax.YGrid = 'off';
c1 = colorbar;%c1.Visible = 'off';
c1.Ticks = hmscale;c1.EdgeColor = 'none';c1.Location = 'westoutside';
c1.TickLabels = {['<=',num2str(hmscale(1))],['>=',num2str(hmscale(2))]};
% OUCH dAIC:motif length = 1:10, realistic g
% [1:568]=-0.540105,[177:276]=-2.039967,[277:376]=6.864094,[377:476]=-2.561496,

%% CDS affinity 2 
fixregion = 1;
% define region relative to ScCrz1
%scregion = [331:430];%NES:186-279;CDS:331-2;CaSc:392;ZF:568
%scregion = [330:429];
binsize = 4;
psite = 0;
motif_length = [1:12];


scregion = cell(1,2);
scregion{1} = [1:568];
%scregion{2} = [229:328]-54;
scregion{2} = [328:427];
%scregion{4} = [429:528]-54;
HMall = nan(2,66);
for i = 1:2
    figure(6),
[CDSA,ID]=cdsaffinity2(folder,fixregion,scregion{i},binsize,psite,motif_length);
HMall(i,:) = CDSA;
Species = flip({ID{1:66}})';
T = table(Species,CDSA);
writetable(T,['Crz1CDSA',num2str(scregion{i}(1)),'to',num2str(scregion{i}(end)),'_v4.csv'],'Delimiter',',','QuoteStrings',true)
end
% heatmap of every region
hmscale = [0 7];
%hmscale = [0 floor(max(HMall(:)))];
figure(7)
HMallsc = HMall; 
HMallsc(HMallsc<hmscale(1)) = hmscale(1); 
HMallsc(HMallsc>hmscale(2)) = hmscale(2);
b=bar3(flip(HMallsc-min(HMallsc(:)),2));%ones(size(Headers)),flip(1:length(Headers)),
for i = 1:66
zdata = b(i).ZData;
b(i).CData = zdata;%(zdata-min(zdata))./(max(zdata)-min(zdata))+1;
b(i).FaceColor = 'interp';b(i).EdgeColor = 'none';
end
view(2)

ax = gca; ax.XColor = 'none';
ax.YColor = 'none';ax.Color = 'none';ax.XGrid = 'off'; ax.YGrid = 'off';
c1 = colorbar;%c1.Visible = 'off';
c1.Ticks = hmscale;c1.EdgeColor = 'none';c1.Location = 'westoutside';
c1.TickLabels = {['<=',num2str(hmscale(1))],['>=',num2str(hmscale(2))]};
% OUCH dAIC:motif length = 1:10, realistic g
% [1:568]=-0.540105,[177:276]=-2.039967,[277:376]=6.864094,[377:476]=-2.561496,


%% CDS affinity change 2 (more Cryptococcus species)
clade = [126,102];%110:Sach clade; 91: Ca clade; 

figure(1)
EdStest = 1;selecON=1;
for W=-3:-1
disp(['W=',num2str(W)])
[S,ance_all,DIS,clade_ID,Z,ID,EdS,nullEdS]=CDSChange2(folder,clade,EdStest,selecON,W);
if selecON ==1, save(['W',num2str(W),'.mat'],'EdS','nullEdS'); end
end
Species = ID(1:66,1);


%% plot CDS across time
Sc = 43;Zr = 55; Kl = 58; Ca = 28; Sp = 7; rest = setdiff(1:66,[Sc,Zr,Kl,Ca,Sp]);
Strains = [Sc,Zr,Kl,Ca,Sp];StX = cell(size(Strains));StY = StX;AAC = StX;
X = [DIS(end)];Y = [ID{end,4}];
for st = 1:length(Strains)
    ance = ance_all{Strains(st)};
    X = [DIS(end)];Y = [ID{end,4}];
    aac = [];ch = [];
    for i = 2:length(ance)

        if isempty(ID{ance(i),2})
            Xs = [X(end), DIS(ance(i))+X(end)];
            Ys = [Y(end), Y(end)];
            X = [X, Xs];
            Y = [Y, Ys];
            aac = [aac, {''}];
        else 
            steps = DIS(ance(i))/(length(ID{ance(i),2})+1);
            changes = ID{ance(i),2};
            aac_all = ID{ance(i),3};
            for j = 1:length(changes)
                Xs = [X(end), steps + X(end)];
                Ys = [changes(j)+Y(end),changes(j)+Y(end)];
                X = [X, Xs];
                Y = [Y, Ys];
                
            end
            aac = [aac, aac_all];
            ch = [ch,changes];
        end
        
    end
    StX{st}=X;
    StY{st}=Y;
    AAC{st}= aac;
    disp(['E[ch]=',num2str(mean(ch))])
end
%
figure(2)
plot(StX{1},StY{1},'r-')
hold on
plot(StX{2},StY{2},'r--')
plot(StX{3},StY{3},'r--')
plot(StX{4},StY{4},'k--')
plot(StX{5},StY{5},'k--')
% for i = 6:38
%     plot(StX{i},StY{i},'k--')
% end
% for i = 39:55
%     plot(StX{i},StY{i},'r-')
% end

% label aa change
for st = 1:length(Strains)
    ance = ance_all{Strains(st)};
    txX = StX{st}; txX = reshape(txX(1:end-1),2,[]);txX = txX(1,:);
    txY = StY{st}; txY = reshape(txY(2:end),2,[]);txY = txY(1,:);
    text(txX,txY+1,AAC{st})
end
hold off
xlabel('Time(AU)')
ylabel({'Predicted calcineurin'; 'docking strength (AU)'})
%% CDS affinity
fixregion = 1;
% define region relative to ScCrz1
scregion = [331:430];%NES:186-279;CDS:331;CaSc:392;ZF:568
binsize = 4;
psite = 0;
figure(3),
[CDSA,ID]=cdsaffinity(folder,fixregion,scregion,binsize,psite);
Species = flip({ID{1:55}})';
%
%CDSA=log(CDSA);
%
T = table(Species,CDSA);
writetable(T,['Crz1CDSA',num2str(scregion(1)),'to',num2str(scregion(end)),'.csv'],'Delimiter',',','QuoteStrings',true)
 %% CDS affinity change 
% clade = [104,85];%104:Sach clade; 85: Ca clade; 105: both
% 
% figure(1)
% [S,ance_all,dS_all,clade_ID,nullEdS,ID]=CDSChange(folder,clade);
% Species = ID(1:55,1);
% %T = table(Species,N,Z);
% %writetable(T,['Crz1CDSAChange.csv'],'Delimiter',',','QuoteStrings',true)
% %%%z test on E[dS] of each clade
% %%%E[dS]Sac=0.25673,p =7.045e-05,z =3.9748
% %%%E[dS]Can=-1.6127,p =0.37832,z =0.881
% %%%E[dS]null=-2.3088+-6.0673 for Sacch,-2.0861+-6.1416 for Ca%% quantify sequences feagure

%% netcharge
fixregion = 1;
% define region relative to ScCrz1
scregion = [392:568];%NES:186-279;CDS:331;CaSc:392;ZF:568
binsize = 3;
figure(2),
netcharge(folder,fixregion,scregion,binsize)
%% netchargedensity
fixregion = 1;
% define region relative to ScCrz1
scregion = [360:568];%NES:186-279;CDS:331;CaSc:392;ZF:568
binsize = 0.05;
psite = 0;
figure(3),
netchargedensity(folder,fixregion,scregion,binsize,psite);
%% sequencechargedecoration
fixregion = 0;
% define region relative to ScCrz1
scregion = [392:568];%NES:186-279;CDS:331;CaSc:392;ZF:568
binsize = 0.15;
psite = 0;
figure(4),
sequencechargedecoration(folder,fixregion,scregion,binsize,psite);
%% sequencechargedecoration sliding window
fixregion = 1;
% define region relative to ScCrz1
scregion = [1:279];%NES:186-279;CDS:331;CaSc:392;ZF:568
binsize = 0.3;
psite = 0;
figure(4),
[~,WR]=sequencechargedecorationSW(folder,fixregion,scregion,binsize,psite);
%%
figure(5),
sequencechargedecorationSW(folder,fixregion,scregion,0.5,1,WR);
%%
figure(3),
netchargedensity(folder,fixregion,scregion,0.025,0,WR);
%% window of lowest SCD
figure(6),
fixregion = 0; 
scregion = [468:568];%NES:186-279;CDS:331;CaSc:392;ZF:568

[SCD0,WR]=sequencechargedecorationSW(folder,fixregion,scregion,0.5,0);
[SCD1,~]=sequencechargedecorationSW(folder,fixregion,scregion,0.5,1,WR);
[NCD0,~]=netchargedensity(folder,fixregion,scregion,0.05,0,WR);
[NCD1,ID]=netchargedensity(folder,fixregion,scregion,0.05,1,WR);
SCDdif = SCD1-SCD0;

Species = flip({ID{1:55}})';
T = table(Species,NCD0,SCD0,SCDdif);
writetable(T,['Crz1SW',num2str(length(scregion)),'.csv'],'Delimiter',',','QuoteStrings',true)
Comparison(NCD0,NCD1,SCD0,SCD1,[1.5,0.05],1,scregion,fixregion,[0,0])
%% fixregion for OUch
figure(6),
fixregion = 1; 
parts = {[1:185],[186:330],[331:568]};%NES:186-279;CDS:331;CaSc:392;ZF:568
for i = 3
    scregion = parts{i};

    [SCD0,~]=sequencechargedecorationSW(folder,fixregion,scregion,0.5,0);
    [SCD1,~]=sequencechargedecorationSW(folder,fixregion,scregion,0.5,1);
    [NCD0,~]=netchargedensity(folder,fixregion,scregion,0.05,0);
    [NCD1,ID]=netchargedensity(folder,fixregion,scregion,0.05,1);
    SCDdif = SCD1-SCD0;
    %
    Species = flip({ID{1:55}})';
    T = table(Species,NCD0,SCD0,SCDdif);
    writetable(T,['Crz1SF',num2str(scregion(1)),'to',num2str(scregion(end)),'.csv'],'Delimiter',',','QuoteStrings',true)
end
%%
Comparison(NCD0,NCD1,SCD0,SCD1,[1.5,0.05],1,scregion,fixregion,[0,0])

%%
function Comparison(NCD0,NCD1,SCD0,SCD1,binsize,Normalize,scregion,fixregion,logON)
ScC = [1:20]; 
CaC = [21:32];
Out = [33:55];%Out = [21:55];

subplot(2,2,1)
plot(NCD0(ScC),SCD0(ScC),'b.',NCD0(CaC),SCD0(CaC),'r.',NCD0(Out),SCD0(Out),'g.')
%plot(NCD0(ScC),SCD0(ScC),'b.',NCD0(Out),SCD0(Out),'r.')
hold on
plot(NCD1(ScC),SCD1(ScC),'bo',NCD1(CaC),SCD1(CaC),'ro',NCD1(Out),SCD1(Out),'go')
%plot(NCD1(ScC),SCD1(ScC),'bo',NCD1(Out),SCD1(Out),'ro')
hold off
text(NCD0(20),SCD0(20),'Sc');text(NCD1(20),SCD1(20),'Sc');
text(NCD0(8),SCD0(8),'Zr');text(NCD1(8),SCD1(8),'Zr');
text(NCD0(6),SCD0(6),'Kl');text(NCD1(6),SCD1(6),'Kl');
text(NCD0(31),SCD0(31),'Ca');text(NCD1(31),SCD1(31),'Ca');
text(NCD0(52),SCD0(52),'Sp');text(NCD1(52),SCD1(52),'Sp');
xlabel('NCD')
ylabel('SCD')
title([num2str(scregion(1)),'-',num2str(scregion(end)),',WithGap=',num2str(fixregion)])

subplot(2,2,2)
NCDdif = NCD1-NCD0;SCDdif = SCD1-SCD0;
if logON(1) == 1, NCDdif = log(NCDdif-min(NCDdif)+1); end
if logON(2) == 1, SCDdif = log(SCDdif-min(SCDdif)+1); end
plot(NCDdif(ScC),SCDdif(ScC),'b.',NCDdif(CaC),SCDdif(CaC),'r.',NCDdif(Out),SCDdif(Out),'g.')
%plot(NCDdif(ScC),SCDdif(ScC),'b.',NCDdif(Out),SCDdif(Out),'r.')
text(NCDdif(20),SCDdif(20),'Sc');
text(NCDdif(8),SCDdif(8),'Zr');
text(NCDdif(6),SCDdif(6),'Kl');
text(NCDdif(31),SCDdif(31),'Ca');
text(NCDdif(52),SCDdif(52),'Sp');
xlabel('Change in NCD')
ylabel('Change in SCD')
save('SCDSWdiff.mat','SCDdif')
subplot(2,2,3)


edges = [min(SCDdif):binsize(1):max(SCDdif)+binsize(1)];
if Normalize ==1
    [NSc]=histcounts(SCDdif(ScC),edges,'Normalization','probability');
    [NCa]=histcounts(SCDdif(CaC),edges,'Normalization','probability');
    [NOut]=histcounts(SCDdif(Out),edges,'Normalization','probability');
    
else
    [NSc]=histcounts(SCDdif(ScC),edges);
    [NCa]=histcounts(SCDdif(CaC),edges);
    [NOut]=histcounts(SCDdif(Out),edges);
    
end
%[NCa]=histcounts(SCDdif(CaC),edges);

bar(edges(2:end),[NSc;NCa;NOut]');%,'stacked')
%bar(edges(2:end),[NSc;NOut]');%,'stacked')
xlabel('Change in SCD')
if Normalize == 1, ylabel('Probability');else,ylabel('Count');end
%ylabel('Count')

subplot(2,2,4)


edges = [min(NCDdif):binsize(2):max(NCDdif)+binsize(2)];
if Normalize ==1
    [NSc]=histcounts(NCDdif(ScC),edges,'Normalization','probability');
    [NCa]=histcounts(NCDdif(CaC),edges,'Normalization','probability');
    [NOut]=histcounts(NCDdif(Out),edges,'Normalization','probability');
    
else
    [NSc]=histcounts(NCDdif(ScC),edges);
    [NCa]=histcounts(NCDdif(CaC),edges);
    [NOut]=histcounts(NCDdif(Out),edges);
    
end
bar(edges(2:end),[NSc;NCa;NOut]');%,'stacked')
%bar(edges(2:end),[NSc;NOut]');%,'stacked')
xlabel('Change in NCD')
if Normalize == 1, ylabel('Probability');else,ylabel('Count');end%ylabel('Count')
end
