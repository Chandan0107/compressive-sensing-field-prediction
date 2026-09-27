close all
clear all
clc

load('Est_mtx_grid_4obj_temp.mat');
load('true_field_4obj_inexact.mat');
load('state_eqn.mat');
load('matrices_data_0.55_25_3.mat');
load('true_tang_fields_dis_5.mat');

% load('meas_10dB_0.3.mat'); % For 10dB;
% snr = 10;

E_true = reshape(E_true,[numel(E_true),1]); E_true(index)=0.01;
E_true = reshape(E_true,[sqrt(length(E_true)),sqrt(length(E_true))]);

method_name = 'nm1';        % nm1 = CS-SOM with data and state constraints (paper Eq. 10)
type        = 'inexact1';
fac = 1; fac1 = 0.7; fac_s = 0.1; M0 = 350;

fname = strcat('error_on_grid_',method_name,'_',num2str(snr),'_0.55','.mat');
delta    = 1e-5;
n_iter = 3;                  % change the n_iter if required
gerror = zeros(n_iter,1);
terror = zeros(n_iter,1);
res_ef = zeros(n_iter,1);
res_sf = zeros(n_iter,1);
xnorm1 = zeros(n_iter,1);
result = [];
result1 = [];
% fprintf('       k  |     p        res_e    res_s    noise_var  res_ef    res_sf    1-norm    t_error   g_error \n')
fprintf('       k  |     res_e    res_s      res_ef    res_sf    t_error   g_error \n')

[U1,S1,V1] = svd(A_state);

for k = 1:n_iter
    b = cell2mat(meas(k));
    A_est = cell2mat(sys_mat(k));
    %     b   = sampled_noise(ff_ob,snr);
    p_n = norm(b)^2/(10^(snr/10));
    tol = sqrt(p_n);
    sig = sqrt(p_n);
    if strcmp(method_name,'nm1')
        %%%%%%%%%%%%%%%%%%%%% NM1 %%%%%%%%%%%%%%%%%%%%
        [U,S,V] = svd(A_est);
        
        %now implement morozov's principle
        %         fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b)/S(p,p);
            res = norm(A_est*xs-b);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        x_svd = xs;
        res_e(k) = norm(A_est*xs - b);
        res_s(k) = norm(A_state*xs - b_state)/norm(b_state);
        %now implement SOM and L1 minimization
        Vs = V(:,1:p);
        %         fac1 = 0.8; fac_s = 0.1;
        Q = eye(2*n)-Vs*Vs'; %the projector matrix onto the other subspace

        is_good_soln = false;
        while ~is_good_soln
            cvx_begin quiet
            variable xn(2*n,1) complex
            minimize(norm(D*(xs+Q*xn),1))
            subject to
            norm(A_est*(xs+Q*xn) - b) <= 1* fac1 * sig;
            norm(A_state*(xs+Q*xn) - b_state)   <= fac_s * norm(b_state);
            cvx_end
            if strcmp(cvx_status,'Solved')
                is_good_soln = true;
            elseif strcmp(cvx_status,'Inaccurate/Solved')
                fac1 = fac1*1.05;
                disp(['Inaccurate/solved, retrying with fac1 = ' num2str(fac1)]);
            else
                %              cvx_stat1{k} = {cvx_status};
                fac1 = fac1*1.05;
                disp(['Failed, retrying with fac1 = ' num2str(fac1)]);
            end
        end
        x_nm= xs + Q*xn;
        
        res_sf(k) = norm(A_state*x_nm-b_state)/norm(b_state);
        
        
    elseif strcmp(method_name,'cs-som')
        %%%%%%%%%%%%%%%%%%%%% CS-SOM %%%%%%%%%%%%%%%%%%%%
        [U,S,V] = svd(A_est);
        
        %now implement morozov's principle
        %         fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b)/S(p,p);
            res = norm(A_est*xs-b);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        x_svd = xs;
        %now implement SOM and L1 minimization
        Vs = V(:,1:p);
        %         fac1 = 0.8;
        Q = eye(2*n)-Vs*Vs'; %the projector matrix onto the other subspace
        is_good_soln = false;
        while ~is_good_soln
            cvx_begin quiet
            variable xn(2*n,1) complex
            minimize(norm(D*(xs+Q*xn),1))
            subject to
            norm(A_est*(xs+Q*xn) - b) <= 1* fac1 * sig;
            cvx_end
            if strcmp(cvx_status,'Solved')
                is_good_soln = true;
            elseif strcmp(cvx_status,'Inaccurate/Solved')
                fac1 = fac1*1.05;
                disp(['Inaccurate/solved, retrying with fac1 = ' num2str(fac1)]);
            else
                %              cvx_stat1{k} = {cvx_status};
                fac1 = fac1*1.05;
                disp(['Failed, retrying with fac1 = ' num2str(fac1)]);
            end
        end
        
        x_nm= xs + Q*xn;
        
    elseif strcmp(method_name,'nm6_norm2')
        
        [U,S,V] = svd(A_est);
        
        %now implement morozov's principle
        %         fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b)/S(p,p);
            res = norm(A_est*xs-b);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        x_svd = xs;
        
        res_e(k) = norm(A_est*xs - b);
        res_s(k) = norm(A_state*xs - b_state)/norm(b_state);
        
        %now implement SOM and L1 minimization
        Vd = V(:,1:p);
        Q = eye(2*n)-Vd*Vd'; %the projector matrix onto the other subspace
        
        Vs = Q*V1(:,1:M0);
        
        %now implement morozov's principle
        %         fac = 1;
        A1 = A_state*Vs;
        b1 = b_state - A_state*xs;
        beta = inv(A1'*A1)*A1'*b1;
        xn   = Vs*beta;
        
        res_ef(k) = norm(A_est*(xs+xn) - b);
        res_sf(k) = norm(A_state*(xs+xn) - b_state)/norm(b_state);
        
        x_nm = xs+xn;
        
        
    elseif strcmp(method_name,'nm6')
        
        [U,S,V] = svd(A_est);
        
        %now implement morozov's principle
        %         fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b)/S(p,p);
            res = norm(A_est*xs-b);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        x_svd = xs;
        
        res_e(k) = norm(A_est*xs - b);
        res_s(k) = norm(A_state*xs - b_state)/norm(b_state);
        
        %now implement SOM and L1 minimization
        Vd = V(:,1:p);
        Q = eye(2*n)-Vd*Vd'; %the projector matrix onto the other subspace
        
        Vs = Q*V1(:,1:M0);
        A1 = A_state*Vs;
        b1 = b_state - A_state*xs;
        beta = inv(A1'*A1)*A1'*b1;
        xn   = Vs*beta;
        %now implement morozov's principle
        %         fac = 1;
        is_good_soln = false;
        while ~is_good_soln
            cvx_begin quiet
            variable xn(M0,1) complex
            minimize(norm(D*(xs+Vs*xn),1))
            subject to
            norm(A_est*(xs+Vs*xn) - b) <= 1* fac1 * sig;
            norm(A_state*(xs+Vs*xn) - b_state)   <= fac_s * norm(b_state);
            cvx_end
            if strcmp(cvx_status,'Solved')
                is_good_soln = true;
            elseif strcmp(cvx_status,'Inaccurate/Solved')
                fac1 = fac1*1.05;
                disp(['Inaccurate/solved, retrying with fac1 = ' num2str(fac1)]);
            else
                %              cvx_stat1{k} = {cvx_status};
                fac1 = fac1*1.05;
                disp(['Failed, retrying with fac1 = ' num2str(fac1)]);
            end
        end
        x_nm= xs + Vs*xn;
        
    elseif strcmp(method_name,'tsvd')
        
        A_n = [A_est];
        b_n = [b];
        
        [U,S,V] = svd(A_n);
        
        %now implement morozov's principle
        %         fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            
            xs = xs + V(:,p) * (U(:,p)'*b_n)/S(p,p);
            res = norm(A_est*xs-b);
            
            if res > fac*sig 
                p = p+1;
            else
                convg = true;
            end
        end
        x_svd = xs;
        
        res_e(k) = norm(A_est*xs - b);
        res_s(k) = norm(A_state*xs - b_state)/norm(b_state);
        
        x_nm = xs;
        
    elseif strcmp(method_name,'nm6_v1')
    
        
        [U,S,V] = svd(A_est);
        
        %now implement morozov's principle
        %         fac = 1;
        p=1; xs = 0; convg = false;
        while ~convg
            xs = xs + V(:,p) * (U(:,p)'*b)/S(p,p);
            res = norm(A_est*xs-b);
            if res > fac*sig
                p = p+1;
            else
                convg = true;
            end
        end
        x_svd = xs;
        
        
        ra = rank(A_state);
        W  = (V(:,1:p)')*V1(:,1:ra);
        [Uw Sw Vw] = svd(W);
        rw = rank(W);
        Wn = Vw(:,rw+1:end);
        
        A1 = A_state*V1(:,1:ra)*Wn;
        b1 = b_state - A_state*xs;
        beta = inv(A1'*A1)*A1'*b1;
        xn   = V1(:,1:ra)*Wn*beta;
        
        
        res_e(k) = norm(A_est*xs - b);
        res_s(k) = norm(A_state*xs - b_state)/norm(b_state);
        
        %now implement SOM and L1 minimization
        x_nm = xs+xn;
   end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    res_ef(k) = norm(A_est*x_nm-b);
    xnorm1(k) = norm(D*x_nm,1);
    terror(k) = norm(x_nm-tang_fields_l)/norm(tang_fields_l);
    res_sf(k) = norm(A_state*x_nm-b_state)/norm(b_state);
    E_pred = B_pred*x_nm;
    E_pred1 = B_pred*x_svd;
    E_pred(index)=0.01;
    E_pred1(index)=0.01;

    E_pred = reshape(E_pred,[sqrt(length(E_pred)),sqrt(length(E_pred))]);
    E_pred1 = reshape(E_pred1,[sqrt(length(E_pred1)),sqrt(length(E_pred1))]);
    gerror(k) = norm(E_pred(3:end-2,3:end-2)-E_true(3:end-2,3:end-2),'fro')/norm(E_true(3:end-2,3:end-2),'fro');
    gerror1(k) = norm(E_pred1(3:end-2,3:end-2)-E_true(3:end-2,3:end-2),'fro')/norm(E_true(3:end-2,3:end-2),'fro');

    % save(fname,'error','Fac1','cvx_stat1');
    result  = [result ;  k   p res_e(k) res_s(k) sig res_ef(k) res_sf(k) xnorm1(k) terror(k) gerror(k) ];
    result1 = [result1 ; k    res_e(k) res_s(k)  res_ef(k) res_sf(k)  terror(k) gerror(k) ];

    disp([result1(k,:)]);
end

mean(result1)

save(fname,'gerror','gerror1','terror','xnorm1','res_sf','res_ef','k','E_pred','E_pred1','E_true','x_nm','tang_fields_l');
% disp(mean(error))
%  disp([' fac,fac1,fac2 '...
%         num2str(fac) ' ' num2str(fac1) ' ' num2str(fac2)  ]);