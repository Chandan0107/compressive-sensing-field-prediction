clear all
close all
clc

load('Closed_Test_4obj.mat');


%% Change the discretization by changing dis
dis=40;

% inexact
Zx=1.25; Zy=1.25; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
[theta1i, rho1i, w1i, N1i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape='rect';
[theta2i, rho2i, w2i, N2i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=1.25; cx=lambda; cy=-2.5*lambda; shape='rect';
[theta3i, rho3i, w3i, N3i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-1*lambda; 
[theta4i, rho4i, w4i, N4i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  shape='rect';
[thetawi, rhowi, wwi, Nwi] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

%% Sample location of each object
Xmw = rhowi.*cos(thetawi); Ymw = rhowi.*sin(thetawi);
Xm1 = rho1i.*cos(theta1i); Ym1 = rho1i.*sin(theta1i);
Xm2 = rho2i.*cos(theta2i); Ym2 = rho2i.*sin(theta2i);
Xm3 = rho3i.*cos(theta3i); Ym3 = rho3i.*sin(theta3i);
Xm4 = rho4i.*cos(theta4i); Ym4 = rho4i.*sin(theta4i);

xx = zeros(N1i,1); yy = zeros(N1i,1);
for i = 1:N1i-1
    xx(i) = (Xm1(i)+Xm1(i+1))/2;
    yy(i) = (Ym1(i)+Ym1(i+1))/2;
end
xx(end) = (Xm1(end)+Xm1(1))/2; Xm1 = xx';
yy(end) = (Ym1(end)+Ym1(1))/2; Ym1 = yy';

xx = zeros(N2i,1); yy = zeros(N2i,1);
for i = 1:N2i-1
    xx(i) = (Xm2(i)+Xm2(i+1))/2;
    yy(i) = (Ym2(i)+Ym2(i+1))/2;
end
xx(end) = (Xm2(end)+Xm2(1))/2; Xm2 = xx';
yy(end) = (Ym2(end)+Ym2(1))/2; Ym2 = yy';

xx = zeros(N3i,1); yy = zeros(N3i,1);
for i = 1:N3i-1
    xx(i) = (Xm3(i)+Xm3(i+1))/2;
    yy(i) = (Ym3(i)+Ym3(i+1))/2;
end
xx(end) = (Xm3(end)+Xm3(1))/2; Xm3 = xx';
yy(end) = (Ym3(end)+Ym3(1))/2; Ym3 = yy';

xx = zeros(N4i,1); yy = zeros(N4i,1);
for i = 1:N4i-1
    xx(i) = (Xm4(i)+Xm4(i+1))/2;
    yy(i) = (Ym4(i)+Ym4(i+1))/2;
end
xx(end) = (Xm4(end)+Xm4(1))/2; Xm4 = xx';
yy(end) = (Ym4(end)+Ym4(1))/2; Ym4 = yy';

xx = zeros(Nwi,1); yy = zeros(Nwi,1);
for i = 1:Nwi-1
    xx(i) = (Xmw(i)+Xmw(i+1))/2;
    yy(i) = (Ymw(i)+Ymw(i+1))/2;
end
xx(end) = (Xmw(end)+Xmw(1))/2; Xmw = xx';
yy(end) = (Ymw(end)+Ymw(1))/2; Ymw = yy';

Xm = [Xm1 Xm2 Xm3 Xm4 Xmw]; 
Ym = [Ym1 Ym2 Ym3 Ym4 Ymw]; 
Nm = length(Xm);

%% Calculate Phi on surfaces 
B_pred = zeros(Nm,2*(N1+N2+N3+N4+Nw)); inc_true=0;
index = [];
for i = 1:Nm
    p = [Xm(i) Ym(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true(i) = incfn(norm(p+0.5*p_ed - Rho));
end

phi_true = B_pred*xcoe1 + inc_true.';



%%%%%%%%%%%%%%%%% gradphi.nhat %%%%%%%%%%%%%%%%%%
del1 = mean(w1i/40);
nx1  = diff([Ym1, Ym1(1)]); nx1 = nx1/max(nx1);
ny1  = -diff([Xm1, Xm1(1)]); ny1 = ny1/max(ny1);

Xm1a = Xm1 - nx1.*del1; Ym1a = Ym1 - ny1.*del1;
Xm1b = Xm1 + nx1.*del1; Ym1b = Ym1 + ny1.*del1;

inc_true1a = 0;
B_pred = zeros(N1i,2*(N1+N2+N3+N4+Nw));
for i = 1:N1i
    p = [Xm1a(i) Ym1a(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true1a(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi1a = B_pred*xcoe1+ inc_true1a.';


inc_true1b = 0;
B_pred = zeros(N1i,2*(N1+N2+N3+N4+Nw));
for i = 1:N1i
    p = [Xm1b(i) Ym1b(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
    
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
    inc_true1b(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi1b = B_pred*xcoe1 + inc_true1b.';
gradphi1_true = ((phi1b-phi1a)/(2*del1));

% 
del2 = mean(w2i/40);
nx2  = diff([Ym2, Ym2(1)]); nx2 = nx2/max(nx2);
ny2  = -diff([Xm2, Xm2(1)]); ny2 = ny2/max(ny2);

Xm2a = Xm2 - nx2.*del2; Ym2a = Ym2 - ny2.*del2;
Xm2b = Xm2 + nx2.*del2; Ym2b = Ym2 + ny2.*del2;

inc_true2a = 0;
B_pred = zeros(N2i,2*(N1+N2+N3+N4+Nw));
for i = 1:N2i
    p = [Xm2a(i) Ym2a(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true2a(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi2a = B_pred*xcoe1 + inc_true2a.';

inc_true2b = 0;
B_pred = zeros(N2i,2*(N1+N2+N3+N4+Nw));
for i = 1:N2i
    p = [Xm2b(i) Ym2b(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true2b(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi2b = B_pred*xcoe1 + inc_true2b.';
gradphi2_true = ((phi2b-phi2a)/(2*del2));

% 
del3 = mean(w3i/40);
nx3  = diff([Ym3, Ym3(1)]); nx3 = nx3/max(nx3);
ny3  = -diff([Xm3, Xm3(1)]); ny3 = ny3/max(ny3);

Xm3a = Xm3 - nx3.*del3; Ym3a = Ym3 - ny3.*del3;
Xm3b = Xm3 + nx3.*del3; Ym3b = Ym3 + ny3.*del3;

inc_true3a = 0;
B_pred = zeros(N3i,2*(N1+N2+N3+N4+Nw));
for i = 1:N3i
    p = [Xm3a(i) Ym3a(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true3a(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi3a = B_pred*xcoe1 + inc_true3a.';

inc_true3b = 0;
B_pred = zeros(N3i,2*(N1+N2+N3+N4+Nw));
for i = 1:N3i
    p = [Xm3b(i) Ym3b(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true3b(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi3b = B_pred*xcoe1 + inc_true3b.';
gradphi3_true = ((phi3b-phi3a)/(2*del3));

%
del4 = mean(w4i/40);
nx4  = diff([Ym4, Ym4(1)]); nx4 = nx4/max(nx4);
ny4  = -diff([Xm4, Xm4(1)]); ny4 = ny4/max(ny4);

Xm4a = Xm4 - nx4.*del4; Ym4a = Ym4 - ny4.*del4;
Xm4b = Xm4 + nx4.*del4; Ym4b = Ym4 + ny4.*del4;

inc_true4a = 0;
B_pred = zeros(N4i,2*(N1+N2+N3+N4+Nw));
for i = 1:N4i
    p = [Xm4a(i) Ym4a(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true4a(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi4a = B_pred*xcoe1 + inc_true4a.';

inc_true4b = 0;
B_pred = zeros(N4i,2*(N1+N2+N3+N4+Nw));
for i = 1:N4i
    p = [Xm4b(i) Ym4b(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_true4b(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phi4b = B_pred*xcoe1+ inc_true4b.';
gradphi4_true = ((phi4b-phi4a)/(2*del4));

%
delw = mean(wwi/40);
nxw  = diff([Ymw, Ymw(1)]); nxw = nxw/max(nxw);
nyw  = -diff([Xmw, Xmw(1)]); nyw = nyw/max(nyw);

Xmwa = Xmw - nxw.*delw; Ymwa = Ymw - nyw.*delw;
Xmwb = Xmw + nxw.*delw; Ymwb = Ymw + nyw.*delw;

inc_truewa = 0;
B_pred = zeros(Nwi,2*(N1+N2+N3+N4+Nw));
for i = 1:Nwi
    p = [Xmwa(i) Ymwa(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_truewa(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phiwa = B_pred*xcoe1 + inc_truewa.';

inc_truewb = 0;
B_pred = zeros(Nwi,2*(N1+N2+N3+N4+Nw));
for i = 1:Nwi
    p = [Xmwb(i) Ymwb(i)];  % ith segment
%     if is_within_contour_pred(p,lambda)
%         index=[index,i];
    
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
    inc_truewb(i) = incfn(norm(p+0.5*p_ed - Rho));
end
phiwb = B_pred*xcoe1 + inc_truewb.';
gradphiw_true = -((phiwb-phiwa)/(2*delw));


% Split the phi fields according to each objects
phi1_true = phi_true(1:N1i);
phi2_true = phi_true(N1i+1 : N1i+N2i);
phi3_true = phi_true(N1i+N2i+1 : N1i+N2i+N3i);
phi4_true = phi_true(N1i+N2i+N3i+1 : N1i+N2i+N3i+N4i);
phiw_true = phi_true(N1i+N2i+N3i+N4i+1 : end);


% save over all fields in sequence 
tang_fields = [gradphi1_true; phi1_true; gradphi2_true;  phi2_true; gradphi3_true; ...
               phi3_true; gradphi4_true;phi4_true;  gradphiw_true; phiw_true] ; 

% plot(abs(tang_fields))
           
save('true_tangfields_all_surfaces_40.mat','tang_fields');