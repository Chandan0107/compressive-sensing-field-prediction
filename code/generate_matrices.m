close all
clear all
clc
tic
load('Closed_Test_4obj.mat');

snr = 25;
dis = 5;                    % discretization lam/5;

% exact
% Zx=1; Zy=1; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
% [theta1i, rho1i, w1i, N1i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx=0.75; Zy=0.75; cx=2*lambda; cy=2*lambda;  shape = 'circ'; 
% [theta2i, rho2i, w2i, N2i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx=2; Zy=0.75; cx=lambda; cy=-2.5*lambda; shape='rect';
% [theta3i, rho3i, w3i, N3i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx=1; Zy=1; cx=-2.5*lambda; cy=-1*lambda; shape = 'circ';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
% [theta4i, rho4i, w4i, N4i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx = 10; Zy = 10; cx = 0; cy = 0;  shape='rect';
% [thetawi, rhowi, wwi, Nwi] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

% inexact
Zx=1.25; Zy=1.25; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
[theta1i, rho1i, w1i, N1i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape='rect';
[theta2i, rho2i, w2i, N2i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=1.25; cx=lambda; cy=-2.5*lambda; shape='rect';
[theta3i, rho3i, w3i, N3i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-1*lambda; 
[theta4i, rho4i, w4i, N4i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

% Zx=1.2; Zy=1.2; cx=-2.5*lambda; cy=-1*lambda; shape='circ'; 
% [theta4i, rho4i, w4i, N4i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  shape='rect';
[thetawi, rhowi, wwi, Nwi] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);


D1=dctmtx(N1i);
D2=dctmtx(N2i);
D3=dctmtx(N3i);
D4=dctmtx(N4i);
Dw=dctmtx(Nwi);
% D2=[];
% D3=[];
% D4=[]; 
% N2i=0; 
% N3i=0; 
% N4i=0;
D=[D1 zeros(N1i,N1i+2*(N2i+N3i+N4i+Nwi));...
    zeros(N1i,N1i) D1 zeros(N1i,2*(N2i+N3i+N4i+Nwi));...
    zeros(N2i,2*(N1i)) D2 zeros(N2i,N2i+2*(N3i+N4i+Nwi));...
    zeros(N2i,2*(N1i)+N2i) D2 zeros(N2i,2*(N3i+N4i+Nwi));...
    zeros(N3i,2*(N1i+N2i)) D3 zeros(N3i,N3i+2*(N4i+Nwi));...
    zeros(N3i,2*(N1i+N2i)+N3i) D3 zeros(N3i,2*(N4i+Nwi));...
    zeros(N4i,2*(N1i+N2i+N3i)) D4 zeros(N4i,N4i+2*(Nwi));...
    zeros(N4i,2*(N1i+N2i+N3i)+N4i) D4 zeros(N4i,2*Nwi);...
    zeros(Nwi,2*(N1i+N2i+N3i+N4i)) Dw zeros(Nwi,Nwi);...
    zeros(Nwi,2*(N1i+N2i+N3i+N4i)+Nwi) Dw];


SR = 0.55;
n    = N1i+N2i+N3i+N4i+Nwi;
Nob=round(SR*2*n); 
n_iter = 3;
meas = cell(n_iter,1);
scat_field = cell(n_iter,1);
sys_mat = cell(n_iter,1);
Px = cell(n_iter,1); Py = cell(n_iter,1);
for k = 1:n_iter
    [P1,P2] = sampling_random(lambda,2*Nob);
    P1=P1(1:Nob); P2=P2(1:Nob);
    ff_ob = zeros(Nob,1);
    Px(k) = {P1}; Py(k) = {P2};
    
    parfor i=1:Nob %observation points
        p = [P1(i) P2(i)];    
        for j=1:N1 %segments of the contour integral
            r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))]; %start of jth segment
            if j~=N1
                r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
            else
                r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            %using one-pt quadrature here as this step is time consumi200ng
            ff_ob(i) =  ff_ob(i) - w1(j) * (xcoe1(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2) ...
                - xcoe1(j+N1)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2));
        end

        for j=1:N2 %segments of the contour integral
            r = [rho2(j)*cos(theta2(j)) rho2(j)*sin(theta2(j))]; %start of jth segment
            if j~=N2
                r_tmp = [rho2(j+1)*cos(theta2(j+1)) rho2(j+1)*sin(theta2(j+1))];
            else
                r_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            %using one-pt quadrature here as this step is time consumi200ng
            ff_ob(i) =  ff_ob(i) - w2(j) * (xcoe1(j+2*(N1))*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2) ...
                - xcoe1(j+2*(N1)+N2)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2));
        end
    %     
        for j=1:N3 %segments of the contour integral
            r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))]; %start of jth segment
            if j~=N3
                r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
            else
                r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            %using one-pt quadrature here as this step is time consumi200ng
            ff_ob(i) =  ff_ob(i) - w3(j) * (xcoe1(j+2*(N1+N2))*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2) ...
                - xcoe1(j+2*(N1+N2)+N3)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2));
        end

        for j=1:N4 %segments of the contour integral
            r = [rho4(j)*cos(theta4(j)) rho4(j)*sin(theta4(j))]; %start of jth segment
            if j~=N4
                r_tmp = [rho4(j+1)*cos(theta4(j+1)) rho4(j+1)*sin(theta4(j+1))];
            else
                r_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            %using one-pt quadrature here as this step is time consumi200ng
            ff_ob(i) =  ff_ob(i) - w4(j) * (xcoe1(j+2*(N1+N2+N3))*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2) ...
                - xcoe1(j+2*(N1+N2+N3)+N4)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2));
        end

        for j=1:Nw %segments of the contour integral
            r = [rhow(j)*cos(thetaw(j)) rhow(j)*sin(thetaw(j))]; %start of jth segment
            if j~=Nw
                r_tmp = [rhow(j+1)*cos(thetaw(j+1)) rhow(j+1)*sin(thetaw(j+1))];
            else
                r_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
            end
            r_ed = r_tmp - r;
            nhat =  -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            %using one-pt quadrature here as this step is time consuming
            ff_ob(i) =  ff_ob(i) - ww(j) * (xcoe1(j+2*(N1+N2+N3+N4))*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2) ...
                - xcoe1(j+2*(N1+N2+N3+N4)+Nw)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2));

        end

    %     ff_ob(i) = ff_ob(i) ; %scattered (0), or total (1) field
        inc(i) = incfn(norm(p-Rho));
    end
    scat_field(k) = {ff_ob};
    b   = sampled_noise(ff_ob,snr);
    meas(k) = {b};
    
    A_est = zeros(Nob,2*n);
    
    for i = 1:Nob
        p = [P1(i) P2(i)];
        % loop over objects
        for j = 1:N1i
            r = [rho1i(j)*cos(theta1i(j)) rho1i(j)*sin(theta1i(j))]; %start of jth segment
            if j~=N1i
                r_tmp = [rho1i(j+1)*cos(theta1i(j+1)) rho1i(j+1)*sin(theta1i(j+1))];
            else
                r_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            A_est(i,j)    = -w1i(j) * glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
            A_est(i,j+N1i) = w1i(j) * glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        end

        for j = 1:N2i
            r = [rho2i(j)*cos(theta2i(j)) rho2i(j)*sin(theta2i(j))]; %start of jth segment
            if j~=N2i
                r_tmp = [rho2i(j+1)*cos(theta2i(j+1)) rho2i(j+1)*sin(theta2i(j+1))];
            else
                r_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            A_est(i,j+2*(N1i))    = -w2i(j) * glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
            A_est(i,j+2*(N1i)+N2i) = w2i(j) * glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        end

        for j = 1:N3i
            r = [rho3i(j)*cos(theta3i(j)) rho3i(j)*sin(theta3i(j))]; %start of jth segment
            if j~=N3i
                r_tmp = [rho3i(j+1)*cos(theta3i(j+1)) rho3i(j+1)*sin(theta3i(j+1))];
            else
                r_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            A_est(i,j+2*(N1i+N2i))    = -w3i(j) * glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
            A_est(i,j+2*(N1i+N2i)+N3i) = w3i(j) * glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        end

        for j = 1:N4i
            r = [rho4i(j)*cos(theta4i(j)) rho4i(j)*sin(theta4i(j))]; %start of jth segment
            if j~=N4i
                r_tmp = [rho4i(j+1)*cos(theta4i(j+1)) rho4i(j+1)*sin(theta4i(j+1))];
            else
                r_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
            end
            r_ed = r_tmp - r;
            nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            A_est(i,j+2*(N1i+N2i+N3i))    = -w4i(j) * glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
            A_est(i,j+2*(N1i+N2i+N3i)+N4i) = w4i(j) * glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        end

        % loop over wall
        for j = 1:Nwi
            r = [rhowi(j)*cos(thetawi(j)) rhowi(j)*sin(thetawi(j))]; %start of jth segment
            if j~=Nwi
                r_tmp = [rhowi(j+1)*cos(thetawi(j+1)) rhowi(j+1)*sin(thetawi(j+1))];
            else
                r_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
            end
            r_ed = r_tmp - r;
            nhat1 =  -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
            A_est(i,j+2*(N1i+N2i+N3i+N4i))     = -wwi(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
            A_est(i,j+2*(N1i+N2i+N3i+N4i)+Nwi)  = wwi(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat1),2);

        end
    end
    sys_mat(k) = {A_est};
    disp(k);
end
fname = strcat('matrices_data_',num2str(SR),'_',num2str(snr),'_',num2str(n_iter),'.mat');
save(fname,'Px','Py','meas','scat_field','sys_mat','n_iter','snr','n','Nob','D');
toc