%-------------------Stress Inversion based on Fractures-------------------%
%                                                                         %
%   program FracInv                                                       %
%   December 2020                                                         %
%                                   By Shihuai Zhang @ ETH Zurich         %
%   Email: zhangshi@ethz.ch                                               %
%          shihuaizhang.xh1@gmail.com                                     %
%-------------------------------------------------------------------------%
clear; close all;
%--------------------------------------------------------------------------
% reading input parameters                                         
%--------------------------------------------------------------------------
run('../Initiation/Input_Fracture')
[depth,dip_dir,dip,criticality]  = textread(input_file,'%f%f%f%f','commentstyle','matlab');
% Fracture Strike
Strike = dip_dir;
for i=1:length(dip_dir)
    if dip_dir(i) < 90
        Strike(i) = dip_dir(i) + 270;
    else
        Strike(i) = dip_dir(i) - 90;
    end
end
%% Prepare for the loop
%--------------------------------------------------------------------------
a_min = 0; a_max = 360; a_int = 45;
b_min = 0; b_max = 90; b_int = 22.5;
c_min = 0; c_max = 90; c_int = 22.5;
phi_min = 0; phi_max = 1; phi_int = 0.2;
S3_min = 0; S3_max = 1; S3_int = 0.5;

a_num = round((a_max-a_min)/a_int);
b_num = round((b_max-b_min)/b_int);
c_num = round((c_max-c_min)/c_int);
phi_num = round((phi_max-phi_min)/phi_int);
S3_num = round((S3_max-S3_min)/S3_int);
total_cnt = a_num*b_num*c_num*phi_num*S3_num;

a_range = zeros(total_cnt,100);
b_range = zeros(total_cnt,100);
c_range = zeros(total_cnt,100);
phi_range = zeros(total_cnt,100);
S3_range = zeros(total_cnt,100);
%--------------------------------------------------------------------------
para_cnt = 0;
for i=1:a_num
    for j=1:b_num
        for m=1:c_num
            for n=1:phi_num
                for h = 1:S3_num
                    para_cnt = para_cnt + 1;
                    a_range(para_cnt,1) = a_min + a_int*(i-1/2);
                    b_range(para_cnt,1) = b_min + b_int*(j-1/2);
                    c_range(para_cnt,1) = c_min + c_int*(m-1/2);
                    phi_range(para_cnt,1) = phi_min + phi_int*(n-1/2);
                    S3_range(para_cnt,1) = S3_min + S3_int*(h-1/2);
                end
            end
        end
    end
end
%% Main loop for the optimal combination of stress parameters
MLE_store = zeros(total_cnt,100); % a, b, c, phi, MLE
iteration_cnt = 0;
while a_int > 0.5
    iteration_cnt = iteration_cnt+1;
    %
    Trial_cnt = 0;
    for lr = 1 : length(a_range)
        Trial_cnt = Trial_cnt + 1;
        a_trial = a_range(lr,iteration_cnt);
        b_trial = b_range(lr,iteration_cnt);
        c_trial = c_range(lr,iteration_cnt);
        phi_trial = phi_range(lr,iteration_cnt);
        S3_trial = S3_range(lr,iteration_cnt);
        % solve stresses on fractures
        [tau, Sn, rake] = Stresses_On_Fractures(Strike,dip,a_trial,b_trial,c_trial,phi_trial,S3_trial);
        % logistic regression
        Slip_Tendency = tau./Sn;
        IS_NaN = isnan(Slip_Tendency(:,1));
        Slip_Tendency(IS_NaN)=[];
        criticality_test = criticality;
        criticality_test(IS_NaN)=[];
        [theta, cost, MLE] = Logistic_Regression(Slip_Tendency,criticality_test);
        % ranking MLE
        MLE_store(Trial_cnt,iteration_cnt) = MLE;
    end
	%
    % select top 40 combinations
    [~,Sort_Index] = sort(MLE_store(:,iteration_cnt),'descend');
    Rank_Top40 = [a_range(Sort_Index(1:40),iteration_cnt)...
        b_range(Sort_Index(1:40),iteration_cnt)...
        c_range(Sort_Index(1:40),iteration_cnt)...
        phi_range(Sort_Index(1:40),iteration_cnt)...
        S3_range(Sort_Index(1:40),iteration_cnt)...
        MLE_store(Sort_Index(1:40),iteration_cnt)];
    %
	% update parameter intervals
    a_int_new = (2/3)*a_int;
    b_int_new = (2/3)*b_int;
    c_int_new = (2/3)*c_int;
    phi_int_new = (2/3)*phi_int;
    S3_int_new = (2/3)*S3_int;
    
    Trial_cnt = 0;
    for ii = 1:2
        for jj = 1:2
            for mm = 1:2
                for nn = 1:2
                    for hh = 1:2
                        for kk = 1:40
                            Trial_cnt = Trial_cnt+1;
                            %
                            a_range(Trial_cnt,iteration_cnt+1) = Rank_Top40(kk,1)-a_int+a_int_new*ii;
                            if a_range(Trial_cnt,iteration_cnt+1) > a_max
                                a_range(Trial_cnt,iteration_cnt+1) = a_max; % parameter should be within the reasonable range
                            end
                            if a_range(Trial_cnt,iteration_cnt+1) < a_min
                                a_range(Trial_cnt,iteration_cnt+1) = a_min;
                            end
                            %      
                            b_range(Trial_cnt,iteration_cnt+1) = Rank_Top40(kk,2)-b_int+b_int_new*jj;
                            if b_range(Trial_cnt,iteration_cnt+1) > b_max
                                b_range(Trial_cnt,iteration_cnt+1) = b_max;
                            end
                            if b_range(Trial_cnt,iteration_cnt+1) < b_min
                                b_range(Trial_cnt,iteration_cnt+1) = b_min;
                            end
                            %
                            c_range(Trial_cnt,iteration_cnt+1) = Rank_Top40(kk,3)-c_int+c_int_new*mm;
                            if abs(b_range(Trial_cnt,iteration_cnt+1)-90)< 5 % assuming Sv is vertical
                                c_range(Trial_cnt,iteration_cnt+1) = 0;
                            end
                            if c_range(Trial_cnt,iteration_cnt+1) > c_max
                                c_range(Trial_cnt,iteration_cnt+1) = c_max;
                            end
                            if c_range(Trial_cnt,iteration_cnt+1) < c_min
                                c_range(Trial_cnt,iteration_cnt+1) = c_min;
                            end
                            %
                            phi_range(Trial_cnt,iteration_cnt+1) = Rank_Top40(kk,4)-phi_int+phi_int_new*nn;
                            if phi_range(Trial_cnt,iteration_cnt+1) > phi_max
                                phi_range(Trial_cnt,iteration_cnt+1) = phi_max;
                            end
                            if phi_range(Trial_cnt,iteration_cnt+1) < phi_min
                                phi_range(Trial_cnt,iteration_cnt+1) = phi_min;
                            end
                            %
                            S3_range(Trial_cnt,iteration_cnt+1) = Rank_Top40(kk,5)-S3_int+S3_int_new*hh;
                            if S3_range(Trial_cnt,iteration_cnt+1) > S3_max
                                S3_range(Trial_cnt,iteration_cnt+1) = S3_max;
                            end
                            if S3_range(Trial_cnt,iteration_cnt+1) < S3_min
                                S3_range(Trial_cnt,iteration_cnt+1) = S3_min;
                            end
                            %
                        end
                    end
                end
            end
        end
    end
    %
    a_int = a_int_new;
    b_int = b_int_new;
    c_int = c_int_new;
    phi_int = phi_int_new;
    S3_int = S3_int_new;
end
%
%
%%
output_file_mat = [output_file,'.mat'];
save(output_file_mat,'MLE_store','a_range','b_range','c_range','phi_range','para_cnt');
