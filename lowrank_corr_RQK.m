function [B,R,RMSE]=lowrank_corr_RQK(X,r,C_i,epsilon,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    GoDec+ Algotithm with RQ kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The code is revised from the code of GoDec+ provided by Xianghao Kong
%INPUTS:
%X: nxp data matrix with n samples and p features
%rank: rank(B)<=rank
%q: >=0, power scheme modification, increasing it lead to better
%OUTPUTS:
%B:BEEG data matrix
%RMSE: Relative error
%R: REEG data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE:
%Xianghao Kong, Wanzeng Kong, Qiaonan Fan, "Task-Independent EEG Identification via Low-Rank Matrix Decomposition"
%Author: Xianghao Kong

%iteration parameters

iter_max=1e+2;
rel_err=[];

X = X';

[m,n]=size(X);


T = zeros(size(X));
B = X;

    iter = 1;
    Y2=randn(n,r);

    while true
        T_sq = T.*T;
        e = T - T.*(1-T_sq./(T_sq+C_i));
        X1=X-e;

            %Update of B
            for i=1:q+1
                Y1=X1*Y2;
                Y2=X1'*Y1;
            end
            [Q,R]=qr(Y2,0);
            base = X1*Q;
            B_new=base*Q';
            Y2 = Q;

        T = X - B_new;

        B_diff = B_new - B;
        stop_cri = (norm(B_diff(:))/norm(B(:)))^2;
        rel_err = [rel_err,stop_cri];

        if stop_cri < epsilon || iter >iter_max
            break;
        end

        B = B_new;
        iter = iter + 1;
    end
    B = B_new;
    R = X-B;
    B = B';
    R = R';
