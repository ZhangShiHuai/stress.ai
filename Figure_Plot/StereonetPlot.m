%%
% Lower-hemisphere stereographic projections of principal stresses
% according to three Euler angles
%
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
%%
[~]=make_stereonet;
hold on
%
% define color
MaxLikelihood = max(max(MLE_store(:,1:Iteration_Num)));
MinLikelihood = min(min(MLE_store(:,1:Iteration_Num)));
rgb1=[0.230, 0.299, 0.754];
rgb2=[0.706, 0.016, 0.150];
C = diverging_map(1000,rgb1,rgb2);
%
for i = 1%:Iteration_Num
    for j=1:para_cnt
        [Sigma_vector_1,Sigma_vector_2,Sigma_vector_3] = EulerAnglesToStressVector(a_range(j,i),b_range(j,i),c_range(j,i));
        %
        Color_ind = round((MLE_store(j,i)-MinLikelihood)/(MaxLikelihood-MinLikelihood)*1000);
        if Color_ind == 0
            Color_ind = 1;
        end
        Color_Array = C(Color_ind,:);
        %
        % Sigma_1
        [x_sigma_1,y_sigma_1] = StressVectorToLowerHemisphere(Sigma_vector_1);
        %
        % Sigma_2
        [x_sigma_2,y_sigma_2] = StressVectorToLowerHemisphere(Sigma_vector_2);
        %
        % Sigma_3
        [x_sigma_3,y_sigma_3] = StressVectorToLowerHemisphere(Sigma_vector_3);
        %
        %
        s1=plot(-y_sigma_1,-x_sigma_1,'>','MarkerEdgeColor',Color_Array,'MarkerSize',10,'LineWidth',1);%i/(i+1)*12
%         maketransparent(s1,round(i/3,1)/17*3)
        hold on
        s2=plot(y_sigma_2,x_sigma_2,'s','MarkerEdgeColor',Color_Array,'MarkerSize',10,'LineWidth',1);
%         maketransparent(s2,round(i/3,1)/17*3)
        hold on
        s3=plot(y_sigma_3,x_sigma_3,'d','MarkerEdgeColor',Color_Array,'MarkerSize',10,'LineWidth',1);
%         maketransparent(s3,round(i/3,1)/17*3)
        hold on
    end
    %
end
%
% set(gcf,'Colormap',C);
% colormap(C);
% caxis([MinLikelihood MaxLikelihood]);
% h = colorbar('FontSize',11,'YTick', linspace(MinLikelihood,MaxLikelihood,5),'YTickLabel',linspace(MinLikelihood,MaxLikelihood,5));
% get(cb,'Position');
% cb.Label.String = 'Likelihood Estimation';
%
%
%
[Sigma_vector_1,Sigma_vector_2,Sigma_vector_3] = EulerAnglesToStressVector(Euler_a,Euler_b,Euler_c);
%
% Sigma_1
[x_sigma_1,y_sigma_1] = StressVectorToLowerHemisphere(Sigma_vector_1);
%
% Sigma_2
[x_sigma_2,y_sigma_2] = StressVectorToLowerHemisphere(Sigma_vector_2);
%
% Sigma_3
[x_sigma_3,y_sigma_3] = StressVectorToLowerHemisphere(Sigma_vector_3);
%
%
plot(-y_sigma_1,-x_sigma_1,'k>','MarkerSize',20,'LineWidth',2.5); % Sigma_1
hold on
plot(y_sigma_2,x_sigma_2,'ks','MarkerSize',20,'LineWidth',2.5); % Sigma_2
hold on
plot(y_sigma_3,x_sigma_3,'kd','MarkerSize',20,'LineWidth',2.5); % Sigma_3
hold on














