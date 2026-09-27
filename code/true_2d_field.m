clear all
close all
clc
% Takes 5 mins
load('Closed_Test_4obj.mat');

lambda=20;
xm   = -4.75*lambda:lambda/20:4.75*lambda;
ym   = fliplr(-4.75*lambda:lambda/20:4.75*lambda);


dis=40;
Zx=1; Zy=1; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
[theta1, rho1, w1, N1] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=0.75; Zy=0.75; cx=2*lambda; cy=2*lambda; shape='circ'; 
[theta2, rho2, w2, N2] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2; Zy=0.75; cx=lambda; cy=-2.5*lambda; shape='rect'; 
[theta3, rho3, w3, N3] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=1; Zy=1; cx=-2.5*lambda; cy=-1*lambda; shape='circ'; 
[theta4, rho4, w4, N4] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx = 10; Zy = 10; cx = 0; cy = 0; shape='rect'; 
[thetaw, rhow, ww, Nw] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

[Xm,Ym] = meshgrid(xm,ym);
Xm = reshape(Xm, [numel(Xm),1]); Ym = reshape(Ym, [numel(Ym),1]);
tic
Nm = length(Xm);
E_true = zeros(Nm,1); E_esti = zeros(Nm,1);
for i = 1:Nm
    p = [Xm(i) Ym(i)];  % ith segment
    A_pred = zeros(1,2*(N1+N2+N3+N4+Nw));
    for j=1:N1 %segments of the contour integral
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))]; %start of jth segment
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat =  [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        A_pred(1,j)    = A_pred(1,j)    - w1(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        A_pred(1,j+N1) = A_pred(1,j+N1) + w1(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        
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
        A_pred(1,j+2*(N1))    = A_pred(1,j+2*(N1))    - w2(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        A_pred(1,j+2*(N1)+N2) = A_pred(1,j+2*(N1)+N2) + w2(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        
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
        A_pred(1,j+2*(N1+N2))    = A_pred(1,j+2*(N1+N2))    - w3(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        A_pred(1,j+2*(N1+N2)+N3) = A_pred(1,j+2*(N1+N2)+N3) + w3(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        
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
        A_pred(1,j+2*(N1+N2+N3))    = A_pred(1,j+2*(N1+N2+N3))    - w4(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        A_pred(1,j+2*(N1+N2+N3)+N4) = A_pred(1,j+2*(N1+N2+N3)+N4) + w4(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
        
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
        A_pred(1,j+2*(N1+N2+N3+N4))    =  A_pred(1,j+2*(N1+N2+N3+N4)) - ww(j)*glquad(@(t)green(k0,p,0,0,r,t,r_ed),2);
        A_pred(1,j+2*(N1+N2+N3+N4)+Nw) =  A_pred(1,j+2*(N1+N2+N3+N4)+Nw) + ww(j)*glquad(@(t)gradgreen(k0,p,0,0,r,t,r_ed,nhat),2);
    end

    E_true(i) =  A_pred*xcoe1; 
    E_esti    = incfn(norm(p-Rho));

end
toc

imagesc(abs(E_true));
figure
% imagesc(abs(E_trues+E_truei));
save('true_field_4obj_inexact.mat','E_true');
