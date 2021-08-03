%% Examplary cell of C.albicans
%cd('C:\Users\Ian\Documents\MATLAB\Crz1\crz1_1_22_2021\Experiment-156\') %0.2 stress
cd('C:\Users\Ian\Documents\MATLAB\Crz1\crz1_2_5_2021\Experiment-159\') % no stress
list = dir(['*t0z0c1*.tif']);
I0 = imread(list.name);
[Ic rect]=imcrop(I0);
imshow(Ic)
%%
n = 3; m = 20;
bigIM=uint8(zeros(n*size(Ic,1),m*size(Ic,2))); 
tp = [1:60/(n*m):60]-1;
k = 0;
for i = 1:n
    for j = 1:m
        k = k+1;
        list = dir(['*t',num2str(tp(k)),'z0c0*.tif']);
        I = imread(list.name);
        Icr = imcrop(I,rect);
        rangerow = (i-1)*size(Icr,1)+1:i*size(Icr,1);
        rangecol = (j-1)*size(Icr,2)+1:j*size(Icr,2);
        bigIM(rangerow,rangecol)=Icr;
        imshow(bigIM)
    end
end
%
imwrite(bigIM,['bigIM',num2str(n),'x',num2str(m),'.tif'],'tif')
%% plot
%cell_id = [37,43];%0.2 stress
cell_id = [58,25,86,92,38];
Ca2M = load('cell_trajectories.mat');
LS_all = Ca2M.LS./Ca2M.PS-1;
figure(2)
plot([ones(length(cell_id),1)*[1:60]/2*60]',LS_all(:,cell_id))
xlim([0 max([1:60]/2*60)])

xlabel('Time (sec)')
ylabel({'Single-cell nuclear';'localization (AU)'})
%% Examplary cell of S.pombe
%cd('C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_8_2021\Prz1_GFP_2021_03_08_3\')%0.2Mstress
cd('C:\Users\Ian\Documents\MATLAB\Crz1\crz1_3_11_2021\Prz1_GFP_2021_03_11_1\')%no stress
list = dir(['*',sprintf('t%03d',1),'z1c2.tif']);
I0 = imread(list.name);
figure(3)
[Ic rect]=imcrop(I0);
imshow(Ic)
%%
n = 2; m = 20;rotate = 1;
bigIM=uint16(zeros(n*size(Ic,1),m*size(Ic,2))); 
if rotate ==1, bigIM=uint16(zeros(n*size(Ic,2),m*size(Ic,1))); end

tp = [301:floor(300/(n*m)):600];
k = 0;
for i = 1:n
    for j = 1:m
        k = k+1;
        list = dir(['*',sprintf('t%03d',tp(k)),'z2c1*.tif']);
        I = imread(list.name);
        Icr = imcrop(I,rect);
        
        if rotate ==1,
            rangerow = (i-1)*size(Icr,2)+1:i*size(Icr,2);
            rangecol = (j-1)*size(Icr,1)+1:j*size(Icr,1);
            bigIM(rangerow,rangecol)=Icr';
        else
            rangerow = (i-1)*size(Icr,1)+1:i*size(Icr,1);
            rangecol = (j-1)*size(Icr,2)+1:j*size(Icr,2);
            bigIM(rangerow,rangecol)=Icr;
        end
            imshow(bigIM)
    end
end
%
imwrite(bigIM,['bigIM',num2str(n),'x',num2str(m),'r',num2str(rotate),'.tif'],'tif')
%%
Pulses = load('PulsesCell#67.mat')
figure(4)
plot([ones(2,1)*[301:601]*6]',Pulses.Pulses.rawLS2Nu(301:601,:))
xlim([300*6 max([1:600]*6)])

xlabel('Time (sec)')
ylabel({'Single-cell nuclear';'localization (AU)'})