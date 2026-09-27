
%% Forward solver use only once. The data is saved as Closed_Test.mat

%% Wall and Object. The object dimension is (lam X lam). 
clear all
close all
clc
tic
%% Set parameters of simulation %
lambda = 20;
k0     = 2*pi/lambda;
tolabs = 1e-4;               % absolute tolerance in integral
tolrel = 1e-2;               % relative tolerance in integral
gamma  = 0.5772; e = exp(1);
eps1    = 3.7 - 2.1i;
eps2    = 1.7 - 1.1i;
eps3    = 2.7 - 3.7i;
eps4    = 1.2 - 0.1i;
epsw    = 3.7 - 2.1i;

k1     = sqrt(eps1)*k0;
k2     = sqrt(eps2)*k0;
k3     = sqrt(eps3)*k0;
k4     = sqrt(eps4)*k0;
kw     = sqrt(epsw)*k0;

%% Object Initialization (square object)
dis = 40;                    % discretization lam/40;

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

%% Scatter
Rho = [lambda/2 -3*lambda/4];
% polarplot(thetaw,rhow); hold on
% polarplot(theta1,rho1); hold on
% polarplot(theta2,rho2); hold on
% polarplot(theta3,rho3); hold on
% polarplot(theta4,rho4); 

%% Matrix Initialization

%N2 = 0;  
% N3 = 0; 
% N4 = 0; 
n   = N1+N2+N3+N4+Nw;
A = zeros(2*n,2*n);
b = zeros(2*n,1);
% The two integrals are:
%oint[g1(r,p) grad(phi).n - grad(g1(r,p).n phi]dr = -phi_inc(p)  (1)
%oint[g2(r,p) grad(phi).n - grad(g2(r,p).n phi]dr = 0            (2)

incfn  = @(rho) besselh(0,2,k0*rho);

gdiag1 = @(k,i) -1j/4*(w1(i) - k^2*w1(i)^3/48 - 1j*(2*w1(i)/pi*(log(w1(i)*k/(4*e))+gamma)));
gdiag2 = @(k,i) -1j/4*(w2(i) - k^2*w2(i)^3/48 - 1j*(2*w2(i)/pi*(log(w2(i)*k/(4*e))+gamma)));
gdiag3 = @(k,i) -1j/4*(w3(i) - k^2*w3(i)^3/48 - 1j*(2*w3(i)/pi*(log(w3(i)*k/(4*e))+gamma)));
gdiag4 = @(k,i) -1j/4*(w4(i) - k^2*w4(i)^3/48 - 1j*(2*w4(i)/pi*(log(w4(i)*k/(4*e))+gamma)));
gdiagw = @(k,i) -1j/4*(ww(i) - k^2*ww(i)^3/48 - 1j*(2*ww(i)/pi*(log(ww(i)*k/(4*e))+gamma)));


for i=1:N1
    p = [rho1(i)*cos(theta1(i)) rho1(i)*sin(theta1(i))];
    if i~=N1
        p_tmp = [rho1(i+1)*cos(theta1(i+1)) rho1(i+1)*sin(theta1(i+1))];
    else
        p_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))];
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A(i,j)       = gdiag1(k0,i);
            A(i,j+N1)    = 1/2.0;
            A(i+N1,j)    = gdiag1(k1,i);
            A(i+N1,j+N1) = -1/2.0;
        else
            A(i,j)       =  w1(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A(i,j+N1)    = -w1(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
            A(i+N1,j)    =  w1(j) * glquad(@(t)green(k1,p,0.5,p_ed,r,t,r_ed),2);
            A(i+N1,j+N1) = -w1(j) * glquad(@(t)gradgreen(k1,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        end
    end
    for j=1:N2
        r = [rho2(j)*cos(theta2(j)) rho2(j)*sin(theta2(j))];
        if j~=N2
            r_tmp = [rho2(j+1)*cos(theta2(j+1)) rho2(j+1)*sin(theta2(j+1))];
        else
            r_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i,j+(2*N1))       =     w2(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i,j+N2+(2*N1))    =    -w2(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        A(i+N1,j+(2*N1))    = 0;
        A(i+N1,j+N2+(2*N1)) = 0;        
    end
    for j=1:N3
        r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))];
        if j~=N3
            r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
        else
            r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i,j+2*(N1+N2))       =     w3(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i,j+N3+2*(N1+N2))    =    -w3(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        A(i+N1,j+2*(N1+N2))    = 0;
        A(i+N1,j+N3+2*(N1+N2)) = 0;       
    end
    for j=1:N4
        r = [rho4(j)*cos(theta4(j)) rho4(j)*sin(theta4(j))];
        if j~=N4
            r_tmp = [rho4(j+1)*cos(theta4(j+1)) rho4(j+1)*sin(theta4(j+1))];
        else
            r_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i,j+2*(N1+N2+N3))       =     w4(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i,j+N4+2*(N1+N2+N3))    =    -w4(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        A(i+N1,j+2*(N1+N2+N3))    = 0;
        A(i+N1,j+N4+2*(N1+N2+N3)) = 0;       
    end
    for j=1:Nw
        r = [rhow(j)*cos(thetaw(j)) rhow(j)*sin(thetaw(j))];
        if j~=Nw
            r_tmp = [rhow(j+1)*cos(thetaw(j+1)) rhow(j+1)*sin(thetaw(j+1))];
        else
            r_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i,j+2*(N1+N2+N3+N4))       =     ww(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i,j+Nw+2*(N1+N2+N3+N4))    =    -ww(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        A(i+N1,j+2*(N1+N2+N3+N4))    = 0;
        A(i+N1,j+Nw+2*(N1+N2+N3+N4)) = 0;
        
    end
    b(i) =  incfn(norm(p+0.5*p_ed - Rho)) ;
end
for i=1:N2
    p = [rho2(i)*cos(theta2(i)) rho2(i)*sin(theta2(i))];
    if i~=N2
        p_tmp = [rho2(i+1)*cos(theta2(i+1)) rho2(i+1)*sin(theta2(i+1))];
    else
        p_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))];
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+(2*N1),j)       =     w1(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+(2*N1),j+N1)    =    -w1(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        A(i+(2*N1)+N2,j)    = 0;
        A(i+(2*N1)+N2,j+N1) = 0;        
    end
    for j=1:N2
        r = [rho2(j)*cos(theta2(j)) rho2(j)*sin(theta2(j))];
        if j~=N2
            r_tmp = [rho2(j+1)*cos(theta2(j+1)) rho2(j+1)*sin(theta2(j+1))];
        else
            r_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A(i+(2*N1),j+(2*N1))       = gdiag2(k0,i);
            A(i+(2*N1),j+(2*N1)+N2)    = 1/2.0;
            A(i+(2*N1)+N2,j+(2*N1))    = gdiag2(k2,i);
            A(i+(2*N1)+N2,j+(2*N1)+N2) = -1/2.0;
        else
            A(i+(2*N1),j+(2*N1))       =  w2(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A(i+(2*N1),j+(2*N1)+N2)    = -w2(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
            A(i+(2*N1)+N2,j+(2*N1))    =  w2(j) * glquad(@(t)green(k2,p,0.5,p_ed,r,t,r_ed),2);
            A(i+(2*N1)+N2,j+(2*N1)+N2) = -w2(j) * glquad(@(t)gradgreen(k2,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        end
    end
    for j=1:N3
        r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))];
        if j~=N3
            r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
        else
            r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+(2*N1),j+2*(N1+N2))       =     w3(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+(2*N1),j+N3+2*(N1+N2))    =    -w3(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        A(i+(2*N1)+N2,j+2*(N1+N2))    = 0;
        A(i+(2*N1)+N2,j+N3+2*(N1+N2)) = 0;       
    end
    for j=1:N4
        r = [rho4(j)*cos(theta4(j)) rho4(j)*sin(theta4(j))];
        if j~=N4
            r_tmp = [rho4(j+1)*cos(theta4(j+1)) rho4(j+1)*sin(theta4(j+1))];
        else
            r_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+(2*N1),j+2*(N1+N2+N3))       =     w4(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+(2*N1),j+N4+2*(N1+N2+N3))    =    -w4(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        A(i+(2*N1)+N2,j+2*(N1+N2+N3))    = 0;
        A(i+(2*N1)+N2,j+N4+2*(N1+N2+N3)) = 0;       
    end
    for j=1:Nw
        r = [rhow(j)*cos(thetaw(j)) rhow(j)*sin(thetaw(j))];
        if j~=Nw
            r_tmp = [rhow(j+1)*cos(thetaw(j+1)) rhow(j+1)*sin(thetaw(j+1))];
        else
            r_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+(2*N1),j+2*(N1+N2+N3+N4))       =     ww(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+(2*N1),j+Nw+2*(N1+N2+N3+N4))    =    -ww(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        A(i+(2*N1)+N2,j+2*(N1+N2+N3+N4))    = 0;
        A(i+(2*N1)+N2,j+Nw+2*(N1+N2+N3+N4)) = 0;
        
    end
    b(i+2*N1) =  incfn(norm(p+0.5*p_ed - Rho)) ;
end
for i=1:N3
    p = [rho3(i)*cos(theta3(i)) rho3(i)*sin(theta3(i))];
    if i~=N3
        p_tmp = [rho3(i+1)*cos(theta3(i+1)) rho3(i+1)*sin(theta3(i+1))];
    else
        p_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))];
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2),j)       =     w1(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2),j+N1)    =    -w1(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        A(i+2*(N1+N2)+N3,j)    = 0;
        A(i+2*(N1+N2)+N3,j+N1) = 0;        
    end
    for j=1:N2
        r = [rho2(j)*cos(theta2(j)) rho2(j)*sin(theta2(j))];
        if j~=N2
            r_tmp = [rho2(j+1)*cos(theta2(j+1)) rho2(j+1)*sin(theta2(j+1))];
        else
            r_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2),j+(2*N1))       =     w2(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2),j+N2+(2*N1))    =    -w2(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        A(i+2*(N1+N2)+N3,j+(2*N1))    = 0;
        A(i+2*(N1+N2)+N3,j+N2+(2*N1)) = 0;        
    end
    for j=1:N3
        r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))];
        if j~=N3
            r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
        else
            r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A(i+2*(N1+N2),j+2*(N1+N2))       = gdiag3(k0,i);
            A(i+2*(N1+N2),j+2*(N1+N2)+N3)    = 1/2.0;
            A(i+2*(N1+N2)+N3,j+2*(N1+N2))    = gdiag3(k3,i);
            A(i+2*(N1+N2)+N3,j+2*(N1+N2)+N3) = -1/2.0;
        else
            A(i+2*(N1+N2),j+2*(N1+N2))       =  w3(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A(i+2*(N1+N2),j+2*(N1+N2)+N3)    = -w3(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
            A(i+2*(N1+N2)+N3,j+2*(N1+N2))    =  w3(j) * glquad(@(t)green(k3,p,0.5,p_ed,r,t,r_ed),2);
            A(i+2*(N1+N2)+N3,j+2*(N1+N2)+N3) = -w3(j) * glquad(@(t)gradgreen(k3,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        end
    end
    for j=1:N4
        r = [rho4(j)*cos(theta4(j)) rho4(j)*sin(theta4(j))];
        if j~=N4
            r_tmp = [rho4(j+1)*cos(theta4(j+1)) rho4(j+1)*sin(theta4(j+1))];
        else
            r_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2),j+2*(N1+N2+N3))       =     w4(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2),j+N4+2*(N1+N2+N3))    =    -w4(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        A(i+2*(N1+N2)+N3,j+2*(N1+N2+N3))    = 0;
        A(i+2*(N1+N2)+N3,j+N4+2*(N1+N2+N3)) = 0;       
    end
    for j=1:Nw
        r = [rhow(j)*cos(thetaw(j)) rhow(j)*sin(thetaw(j))];
        if j~=Nw
            r_tmp = [rhow(j+1)*cos(thetaw(j+1)) rhow(j+1)*sin(thetaw(j+1))];
        else
            r_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2),j+2*(N1+N2+N3+N4))       =     ww(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2),j+Nw+2*(N1+N2+N3+N4))    =    -ww(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        A(i+2*(N1+N2)+N3,j+2*(N1+N2+N3+N4))    = 0;
        A(i+2*(N1+N2)+N3,j+Nw+2*(N1+N2+N3+N4)) = 0;
        
    end
    b(i+2*(N1+N2)) =  incfn(norm(p+0.5*p_ed - Rho)) ;
end
for i=1:N4
    p = [rho4(i)*cos(theta4(i)) rho4(i)*sin(theta4(i))];
    if i~=N4
        p_tmp = [rho4(i+1)*cos(theta4(i+1)) rho4(i+1)*sin(theta4(i+1))];
    else
        p_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))];
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3),j)       =     w1(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3),j+N1)    =    -w1(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        A(i+2*(N1+N2+N3)+N4,j)    = 0;
        A(i+2*(N1+N2+N3)+N4,j+N1) = 0;        
    end
    for j=1:N2
        r = [rho2(j)*cos(theta2(j)) rho2(j)*sin(theta2(j))];
        if j~=N2
            r_tmp = [rho2(j+1)*cos(theta2(j+1)) rho2(j+1)*sin(theta2(j+1))];
        else
            r_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3),j+(2*N1))       =     w2(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3),j+N2+(2*N1))    =    -w2(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        A(i+2*(N1+N2+N3)+N4,j+(2*N1))    = 0;
        A(i+2*(N1+N2+N3)+N4,j+N2+(2*N1)) = 0;        
    end
    for j=1:N3
        r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))];
        if j~=N3
            r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
        else
            r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3),j+2*(N1+N2))       =     w3(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3),j+N3+2*(N1+N2))    =    -w3(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        A(i+2*(N1+N2+N3)+N4,j+2*(N1+N2))    = 0;
        A(i+2*(N1+N2+N3)+N4,j+N3+2*(N1+N2)) = 0;       
    end
    for j=1:N4
        r = [rho4(j)*cos(theta4(j)) rho4(j)*sin(theta4(j))];
        if j~=N4
            r_tmp = [rho4(j+1)*cos(theta4(j+1)) rho4(j+1)*sin(theta4(j+1))];
        else
            r_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A(i+2*(N1+N2+N3),j+2*(N1+N2+N3))       = gdiag4(k0,i);
            A(i+2*(N1+N2+N3),j+2*(N1+N2+N3)+N4)    = 1/2.0;
            A(i+2*(N1+N2+N3)+N4,j+2*(N1+N2+N3))    = gdiag4(k4,i);
            A(i+2*(N1+N2+N3)+N4,j+2*(N1+N2+N3)+N4) = -1/2.0;
        else
            A(i+2*(N1+N2+N3),j+2*(N1+N2+N3))       =  w4(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A(i+2*(N1+N2+N3),j+2*(N1+N2+N3)+N4)    = -w4(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
            A(i+2*(N1+N2+N3)+N4,j+2*(N1+N2+N3))    =  w4(j) * glquad(@(t)green(k4,p,0.5,p_ed,r,t,r_ed),2);
            A(i+2*(N1+N2+N3)+N4,j+2*(N1+N2+N3)+N4) = -w4(j) * glquad(@(t)gradgreen(k4,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        end
    end
    for j=1:Nw
        r = [rhow(j)*cos(thetaw(j)) rhow(j)*sin(thetaw(j))];
        if j~=Nw
            r_tmp = [rhow(j+1)*cos(thetaw(j+1)) rhow(j+1)*sin(thetaw(j+1))];
        else
            r_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3),j+2*(N1+N2+N3+N4))       =     ww(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3),j+Nw+2*(N1+N2+N3+N4))    =    -ww(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        A(i+2*(N1+N2+N3)+N4,j+2*(N1+N2+N3+N4))    = 0;
        A(i+2*(N1+N2+N3)+N4,j+Nw+2*(N1+N2+N3+N4)) = 0;
        
    end
    b(i+2*(N1+N2+N3)) =  incfn(norm(p+0.5*p_ed - Rho)) ;
end
for i=1:Nw
    p = [rhow(i)*cos(thetaw(i)) rhow(i)*sin(thetaw(i))];
    if i~=Nw
        p_tmp = [rhow(i+1)*cos(thetaw(i+1)) rhow(i+1)*sin(thetaw(i+1))];
    else
        p_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1
        r = [rho1(j)*cos(theta1(j)) rho1(j)*sin(theta1(j))];
        if j~=N1
            r_tmp = [rho1(j+1)*cos(theta1(j+1)) rho1(j+1)*sin(theta1(j+1))];
        else
            r_tmp = [rho1(1)*cos(theta1(1)) rho1(1)*sin(theta1(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3+N4),j)       =     w1(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3+N4),j+N1)    =    -w1(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        A(i+2*(N1+N2+N3+N4)+Nw,j)    = 0;
        A(i+2*(N1+N2+N3+N4)+Nw,j+N1) = 0;        
    end
    for j=1:N2
        r = [rho2(j)*cos(theta2(j)) rho2(j)*sin(theta2(j))];
        if j~=N2
            r_tmp = [rho2(j+1)*cos(theta2(j+1)) rho2(j+1)*sin(theta2(j+1))];
        else
            r_tmp = [rho2(1)*cos(theta2(1)) rho2(1)*sin(theta2(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3+N4),j+(2*N1))       =     w2(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3+N4),j+N2+(2*N1))    =    -w2(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        A(i+2*(N1+N2+N3+N4)+Nw,j+(2*N1))    = 0;
        A(i+2*(N1+N2+N3+N4)+Nw,j+N2+(2*N1)) = 0;        
    end
    for j=1:N3
        r = [rho3(j)*cos(theta3(j)) rho3(j)*sin(theta3(j))];
        if j~=N3
            r_tmp = [rho3(j+1)*cos(theta3(j+1)) rho3(j+1)*sin(theta3(j+1))];
        else
            r_tmp = [rho3(1)*cos(theta3(1)) rho3(1)*sin(theta3(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3+N4),j+2*(N1+N2))       =     w3(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3+N4),j+N3+2*(N1+N2))    =    -w3(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        A(i+2*(N1+N2+N3+N4)+Nw,j+2*(N1+N2))    = 0;
        A(i+2*(N1+N2+N3+N4)+Nw,j+N3+2*(N1+N2)) = 0;       
    end
    for j=1:N4
        r = [rho4(j)*cos(theta4(j)) rho4(j)*sin(theta4(j))];
        if j~=N4
            r_tmp = [rho4(j+1)*cos(theta4(j+1)) rho4(j+1)*sin(theta4(j+1))];
        else
            r_tmp = [rho4(1)*cos(theta4(1)) rho4(1)*sin(theta4(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A(i+2*(N1+N2+N3+N4),j+2*(N1+N2+N3))       =     w4(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A(i+2*(N1+N2+N3+N4),j+N4+2*(N1+N2+N3))    =    -w4(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        A(i+2*(N1+N2+N3+N4)+Nw,j+2*(N1+N2+N3))    = 0;
        A(i+2*(N1+N2+N3+N4)+Nw,j+N4+2*(N1+N2+N3)) = 0;       
    end
    for j=1:Nw
        r = [rhow(j)*cos(thetaw(j)) rhow(j)*sin(thetaw(j))];
        if j~=Nw
            r_tmp = [rhow(j+1)*cos(thetaw(j+1)) rhow(j+1)*sin(thetaw(j+1))];
        else
            r_tmp = [rhow(1)*cos(thetaw(1)) rhow(1)*sin(thetaw(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A(i+2*(N1+N2+N3+N4),j+2*(N1+N2+N3+N4))       = gdiagw(k0,i);
            A(i+2*(N1+N2+N3+N4),j+2*(N1+N2+N3+N4)+Nw)    = 1/2.0;
            A(i+2*(N1+N2+N3+N4)+Nw,j+2*(N1+N2+N3+N4))    = gdiagw(kw,i);
            A(i+2*(N1+N2+N3+N4)+Nw,j+2*(N1+N2+N3+N4)+Nw) = -1/2.0;
        else
            A(i+2*(N1+N2+N3+N4),j+2*(N1+N2+N3+N4))       =  ww(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A(i+2*(N1+N2+N3+N4),j+2*(N1+N2+N3+N4)+Nw)    = -ww(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
            A(i+2*(N1+N2+N3+N4)+Nw,j+2*(N1+N2+N3+N4))    =  ww(j) * glquad(@(t)green(kw,p,0.5,p_ed,r,t,r_ed),2);
            A(i+2*(N1+N2+N3+N4)+Nw,j+2*(N1+N2+N3+N4)+Nw) = -ww(j) * glquad(@(t)gradgreen(kw,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        end
    end
    b(i+2*(N1+N2+N3+N4)) =  incfn(norm(p+0.5*p_ed-Rho)) ;
end

% Solve the matrix
xcoe1 = A\b;
% figure
% plot(abs(xcoe1))
save('Closed_Test_4obj.mat');
toc