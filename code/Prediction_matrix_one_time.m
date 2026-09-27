clear all
close all
clc
% Takes 4 mins
load('Closed_Test_4obj.mat','lambda','k0'); 
dis = 5;                    % discretization lam/5;

% exact
% Zx=1; Zy=1; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
% [theta1, rho1, w1, N1] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx=0.75; Zy=0.75; cx=2*lambda; cy=2*lambda;  shape = 'circ'; 
% [theta2, rho2, w2, N2] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx=2; Zy=0.75; cx=lambda; cy=-2.5*lambda;  shape='rect';
% [theta3, rho3, w3, N3] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx=1; Zy=1; cx=-2.5*lambda; cy=-1*lambda; shape = 'circ';
% [theta4, rho4, w4, N4] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx = 10; Zy = 10; cx = 0; cy = 0; shape='rect'; 
% [thetaw, rhow, ww, Nw] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

% inexact
Zx=1.25; Zy=1.25; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
[theta1, rho1, w1, N1] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape = 'rect';
[theta2, rho2, w2, N2] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=1.25; cx=lambda; cy=-2.5*lambda; shape = 'rect';  
[theta3, rho3, w3, N3] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

% Zx=2.2; Zy=2.2; cx=-2.5*lambda; cy=-1*lambda; 
% [theta4, rho4, w4, N4] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-1*lambda; 
[theta4, rho4, w4, N4] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  shape = 'rect';
[thetaw, rhow, ww, Nw] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);



% N2 = 0; 
% N3 = 0; 
% N4 = 0;

xm   = -4.75*lambda:lambda/20:4.75*lambda;
ym   = fliplr(-4.75*lambda:lambda/20:4.75*lambda);
[Xm,Ym] = meshgrid(xm,ym);
Xm = reshape(Xm, [numel(Xm),1]); Ym = reshape(Ym, [numel(Ym),1]);
% obs = 2.5*lambda;
% obs_n = 301; obs_t = 2*pi/(obs_n-1);
% rxth = 0:obs_t:2*pi;
% Xm   = [obs.*cos(rxth)];
% Ym   = [obs.*sin(rxth)];
Nm = length(Xm);
tic
B_pred = zeros(Nm,2*(N1+N2+N3+N4+Nw));
index = [];
for i = 1:Nm
    p = [Xm(i) Ym(i)];  % ith segment
    if is_within_contour_pred(p,lambda)
        index=[index,i];
    end
    for j=1:N1 %segments of the contour integral
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))]; %start of jth segment
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        B_pred(i,j)    = B_pred(i,j)    - w1(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        B_pred(i,j+N1) = B_pred(i,j+N1) + w1(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);

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
        B_pred(i,j+2*(N1))    = B_pred(i,j+2*(N1))    - w2(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        B_pred(i,j+2*(N1)+N2) = B_pred(i,j+2*(N1)+N2) + w2(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);

    end
    
    for j=1:N3 %segments of the contour integral
        r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))]; %start of jth segment
        if j~=N3
            r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
        else
            r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
        end
        r_ed = r_tmp - r;
        nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        B_pred(i,j+2*(N1+N2))    = B_pred(i,j+2*(N1+N2))    - w3(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        B_pred(i,j+2*(N1+N2)+N3) = B_pred(i,j+2*(N1+N2)+N3) + w3(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);

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
        B_pred(i,j+2*(N1+N2+N3))    = B_pred(i,j+2*(N1+N2+N3))    - w4(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        B_pred(i,j+2*(N1+N2+N3)+N4) = B_pred(i,j+2*(N1+N2+N3)+N4) + w4(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);

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
        B_pred(i,j+2*(N1+N2+N3+N4))    =  B_pred(i,j+2*(N1+N2+N3+N4)) - ww(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        B_pred(i,j+2*(N1+N2+N3+N4)+Nw) =  B_pred(i,j+2*(N1+N2+N3+N4)+Nw) + ww(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
    end
end
toc
save('Est_mtx_grid_4obj_temp.mat','B_pred','index');
beep