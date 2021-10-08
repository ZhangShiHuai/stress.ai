%% Plot iterations of each stress parameter
clear
run('../Initiation/Input_Fracture')
load([output_file '.mat'])
%% Field measurements
Euler_a = 237; % or 237 - 180 = 57
Euler_b = 0;
Euler_c = 90;
Phi = 0.6022;
%% delete zero-columns
a_range(:,all(a_range==0,1)) = [];
b_range(:,all(b_range==0,1)) = [];
c_range(:,all(c_range==0,1)) = [];
phi_range(:,all(phi_range==0,1)) = [];
MLE_store(:,all(MLE_store==0,1)) = [];
%
[~,Iteration_Num] = size(MLE_store);
%% Euler angle a
% define color
MaxLikelihood = max(max(MLE_store(:,1:Iteration_Num)));
MinLikelihood = min(min(MLE_store(:,1:Iteration_Num)));
rgb1=[0.230, 0.299, 0.754];
rgb2=[0.706, 0.016, 0.150];
C = diverging_map(256,rgb1,rgb2);
%%
figure
hold on
plot([0 Iteration_Num],[Euler_a Euler_a],'r:','LineWidth',2) % Cajon Pass
%
xlim([0,Iteration_Num])
ylim([0,360])
set(gca,'ytick',0:45:360)
set(gca,'xtick',1:1:Iteration_Num)
grid on
box on
%
%
for i = 1:Iteration_Num
    for j=1:(para_cnt-1)
        if a_range(j,i)==a_range(j+1,i)
        else
            Color_ind = round((MLE_store(j,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
            Color_Array = C(Color_ind,:);
            plot(i,a_range(j,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
            hold on
        end
    end
    Color_ind = round((MLE_store(j+1,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
    Color_Array = C(Color_ind,:);
    plot(i,a_range(j+1,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
    hold on
end
%
set(gcf,'Colormap',C);
colormap(C);
caxis([MinLikelihood MaxLikelihood]);
cb = colorbar;
get(cb,'Position');
cb.Label.String = 'Likelihood Estimation';
%
%% Euler angle b
figure
plot([0 Iteration_Num],[Euler_b Euler_b],'r:','LineWidth',2) % Cajon Pass
%
xlim([0,Iteration_Num])
ylim([0,90])
set(gca,'ytick',0:30:90)
set(gca,'xtick',1:1:Iteration_Num)
grid on
box on
%
hold on
for i = 1:Iteration_Num
    for j=1:(para_cnt-1)
        if b_range(j,i)==b_range(j+1,i)
        else
            Color_ind = round((MLE_store(j,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
            Color_Array = C(Color_ind,:);
            plot(i,b_range(j,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
            hold on
        end
    end
    Color_ind = round((MLE_store(j+1,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
    Color_Array = C(Color_ind,:);
    plot(i,b_range(j+1,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
    hold on
end
hold on

%% Euler angle c
figure
plot([0 Iteration_Num],[Euler_c Euler_c],'r:','LineWidth',2) % Cajon Pass
%
xlim([0,Iteration_Num])
ylim([0,90])
set(gca,'ytick',0:30:90)
set(gca,'xtick',1:1:Iteration_Num)
grid on
box on
%
hold on
for i = 1:Iteration_Num
    for j=1:(para_cnt-1)
        if c_range(j,i)==c_range(j+1,i)
        else
            Color_ind = round((MLE_store(j,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
            Color_Array = C(Color_ind,:);
            plot(i,c_range(j,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
            hold on
        end
    end
    Color_ind = round((MLE_store(j+1,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
    Color_Array = C(Color_ind,:);
    plot(i,c_range(j+1,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
    hold on
end
hold on

%%  Euler angle phi
figure
plot([0 Iteration_Num],[Phi Phi],'r:','LineWidth',2) % Cajon Pass
%
xlim([0,Iteration_Num])
ylim([0,1])
set(gca,'ytick',0:0.1:1)
set(gca,'xtick',1:1:Iteration_Num)
grid on
box on
hold on
%
for i = 1:Iteration_Num
    for j=1:(para_cnt-1)
        if phi_range(j,i)==phi_range(j+1,i)
        else
            Color_ind = round((MLE_store(j,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
            Color_Array = C(Color_ind,:);
            plot(i,phi_range(j,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
            hold on
        end
    end
    Color_ind = round((MLE_store(j+1,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*255)+1;
    Color_Array = C(Color_ind,:);
    plot(i,phi_range(j+1,i),'o','Color',Color_Array,'MarkerFaceColor',Color_Array)
    hold on
end
hold on
