close all
clear all
clc
load('Closed_Test_4obj.mat');

dis = 5;                     % discretization lam/5;

% %exact
% Zx=1; Zy=1; cx=-2*lambda; cy=1.5*lambda; shape='rect'; 
% [theta1i, rho1i, w1i, N1i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);
% 
% Zx = 10; Zy = 10; cx = 0; cy = 0; shape='rect'; 
% [thetawi, rhowi, wwi, Nwi] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

% inexact
Zx = 1.25; Zy = 1.25; cx=-2*lambda; cy=2.5*lambda; shape='rect'; 
[theta1i, rho1i, w1i, N1i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape='rect';
[theta2i, rho2i, w2i, N2i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=1.25; cx=2*lambda; cy=-2.5*lambda; shape='rect';
[theta3i, rho3i, w3i, N3i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-2*lambda; 
[theta4i, rho4i, w4i, N4i] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  shape='rect';
[thetawi, rhowi, wwi, Nwi] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);


% diagonal element for state equation
gdiag1s = @(k,i) 1j/4*(w1i(i) - k^2*w1i(i)^3/48 - 1j*(2*w1i(i)/pi*(log(w1i(i)*k/(4*e))+gamma)));
gdiag2s = @(k,i) 1j/4*(w2i(i) - k^2*w2i(i)^3/48 - 1j*(2*w2i(i)/pi*(log(w2i(i)*k/(4*e))+gamma)));
gdiag3s = @(k,i) 1j/4*(w3i(i) - k^2*w3i(i)^3/48 - 1j*(2*w3i(i)/pi*(log(w3i(i)*k/(4*e))+gamma)));
gdiag4s = @(k,i) 1j/4*(w4i(i) - k^2*w4i(i)^3/48 - 1j*(2*w4i(i)/pi*(log(w4i(i)*k/(4*e))+gamma)));
gdiagws = @(k,i) 1j/4*(wwi(i) - k^2*wwi(i)^3/48 - 1j*(2*wwi(i)/pi*(log(wwi(i)*k/(4*e))+gamma)));


% N2i = 0;
% N3i = 0;
% N4i = 0;

n    = N1i+N2i+N3i+N4i+Nwi;

% Populating state eqn coefficients
A_state = zeros(n,2*n);
b_state = zeros(n,1);

for i=1:N1i
    p = [rho1i(i)*cos(theta1i(i)) rho1i(i)*sin(theta1i(i))];
    if i~=N1i
        p_tmp = [rho1i(i+1)*cos(theta1i(i+1)) rho1i(i+1)*sin(theta1i(i+1))];
    else
        p_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1i
        r = [rho1i(j)*cos(theta1i(j)) rho1i(j)*sin(theta1i(j))];
        if j~=N1i
            r_tmp = [rho1i(j+1)*cos(theta1i(j+1)) rho1i(j+1)*sin(theta1i(j+1))];
        else
            r_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A_state(i,j)       = gdiag1s(k0,i);
            A_state(i,j+N1i)   = -1/2.0; % It is -1/2 as a_i is brought to LHS
        else
            A_state(i,j)       =  -w1i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A_state(i,j+N1i)   = w1i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        end
    end
    for j=1:N2i
        r = [rho2i(j)*cos(theta2i(j)) rho2i(j)*sin(theta2i(j))];
        if j~=N2i
            r_tmp = [rho2i(j+1)*cos(theta2i(j+1)) rho2i(j+1)*sin(theta2i(j+1))];
        else
            r_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i,j+(2*N1i))       =     -w2i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i,j+N2i+(2*N1i))    =    w2i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);        
    end   
    for j=1:N3i
        r = [rho3i(j)*cos(theta3i(j)) rho3i(j)*sin(theta3i(j))];
        if j~=N3i
            r_tmp = [rho3i(j+1)*cos(theta3i(j+1)) rho3i(j+1)*sin(theta3i(j+1))];
        else
            r_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i,j+2*(N1i+N2i))       =     -w3i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i,j+N3i+2*(N1i+N2i))    =    w3i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);       
    end
    for j=1:N4i
        r = [rho4i(j)*cos(theta4i(j)) rho4i(j)*sin(theta4i(j))];
        if j~=N4i
            r_tmp = [rho4i(j+1)*cos(theta4i(j+1)) rho4i(j+1)*sin(theta4i(j+1))];
        else
            r_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i,j+2*(N1i+N2i+N3i))       =     -w4i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i,j+N4i+2*(N1i+N2i+N3i))    =    w4i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
    end
    for j=1:Nwi
        r = [rhowi(j)*cos(thetawi(j)) rhowi(j)*sin(thetawi(j))];
        if j~=Nwi
            r_tmp = [rhowi(j+1)*cos(thetawi(j+1)) rhowi(j+1)*sin(thetawi(j+1))];
        else
            r_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i,j+2*(N1i+N2i+N3i+N4i))        =     -wwi(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i,j+Nwi+2*(N1i+N2i+N3i+N4i))    =    wwi(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
    end
    
    b_state(i) =  -incfn(norm(p+0.5*p_ed - Rho)) ;
end
for i=1:N2i
    p = [rho2i(i)*cos(theta2i(i)) rho2i(i)*sin(theta2i(i))];
    if i~=N2i
        p_tmp = [rho2i(i+1)*cos(theta2i(i+1)) rho2i(i+1)*sin(theta2i(i+1))];
    else
        p_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1i
        r = [rho1i(j)*cos(theta1i(j)) rho1i(j)*sin(theta1i(j))];
        if j~=N1i
            r_tmp = [rho1i(j+1)*cos(theta1i(j+1)) rho1i(j+1)*sin(theta1i(j+1))];
        else
            r_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+N1i,j)       =     -w1i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+N1i,j+N1i)    =    w1i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
    end
    for j=1:N2i
        r = [rho2i(j)*cos(theta2i(j)) rho2i(j)*sin(theta2i(j))];
        if j~=N2i
            r_tmp = [rho2i(j+1)*cos(theta2i(j+1)) rho2i(j+1)*sin(theta2i(j+1))];
        else
            r_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A_state(i+(N1i),j+(2*N1i))       = gdiag2s(k0,i);
            A_state(i+(N1i),j+(2*N1i)+N2i)    = -1/2.0;
        else
            A_state(i+N1i,j+(2*N1i))       =  -w2i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A_state(i+N1i,j+(2*N1i)+N2i)    = w2i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        end
    end
    for j=1:N3i
        r = [rho3i(j)*cos(theta3i(j)) rho3i(j)*sin(theta3i(j))];
        if j~=N3i
            r_tmp = [rho3i(j+1)*cos(theta3i(j+1)) rho3i(j+1)*sin(theta3i(j+1))];
        else
            r_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+N1i,j+2*(N1i+N2i))       =     -w3i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+N1i,j+N3i+2*(N1i+N2i))    =    w3i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
    end
    for j=1:N4i
        r = [rho4i(j)*cos(theta4i(j)) rho4i(j)*sin(theta4i(j))];
        if j~=N4i
            r_tmp = [rho4i(j+1)*cos(theta4i(j+1)) rho4i(j+1)*sin(theta4i(j+1))];
        else
            r_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i),j+2*(N1i+N2i+N3i))       =     -w4i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i),j+N4i+2*(N1i+N2i+N3i))    =    w4i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
    end
    for j=1:Nwi
        r = [rhowi(j)*cos(thetawi(j)) rhowi(j)*sin(thetawi(j))];
        if j~=Nwi
            r_tmp = [rhowi(j+1)*cos(thetawi(j+1)) rhowi(j+1)*sin(thetawi(j+1))];
        else
            r_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i),j+2*(N1i+N2i+N3i+N4i))       =     -wwi(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i),j+Nwi+2*(N1i+N2i+N3i+N4i))    =    wwi(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
       
    end
    b_state(i+N1i) =  -incfn(norm(p+0.5*p_ed - Rho)) ;

end
for i=1:N3i
    p = [rho3i(i)*cos(theta3i(i)) rho3i(i)*sin(theta3i(i))];
    if i~=N3i
        p_tmp = [rho3i(i+1)*cos(theta3i(i+1)) rho3i(i+1)*sin(theta3i(i+1))];
    else
        p_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1i
        r = [rho1i(j)*cos(theta1i(j)) rho1i(j)*sin(theta1i(j))];
        if j~=N1i
            r_tmp = [rho1i(j+1)*cos(theta1i(j+1)) rho1i(j+1)*sin(theta1i(j+1))];
        else
            r_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i),j)       =     -w1i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i),j+N1i)    =    w1i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
    end
    for j=1:N2i
        r = [rho2i(j)*cos(theta2i(j)) rho2i(j)*sin(theta2i(j))];
        if j~=N2i
            r_tmp = [rho2i(j+1)*cos(theta2i(j+1)) rho2i(j+1)*sin(theta2i(j+1))];
        else
            r_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i),j+(2*N1i))       =     -w2i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i),j+N2i+(2*N1i))    =    w2i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
       
    end
    for j=1:N3i
        r = [rho3i(j)*cos(theta3i(j)) rho3i(j)*sin(theta3i(j))];
        if j~=N3i
            r_tmp = [rho3i(j+1)*cos(theta3i(j+1)) rho3i(j+1)*sin(theta3i(j+1))];
        else
            r_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A_state(i+(N1i+N2i),j+2*(N1i+N2i))       = gdiag3s(k0,i);
            A_state(i+(N1i+N2i),j+2*(N1i+N2i)+N3i)    = -1/2.0;
            
        else
            A_state(i+(N1i+N2i),j+2*(N1i+N2i))       =  -w3i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A_state(i+(N1i+N2i),j+2*(N1i+N2i)+N3i)    = w3i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
            
        end
    end
    for j=1:N4i
        r = [rho4i(j)*cos(theta4i(j)) rho4i(j)*sin(theta4i(j))];
        if j~=N4i
            r_tmp = [rho4i(j+1)*cos(theta4i(j+1)) rho4i(j+1)*sin(theta4i(j+1))];
        else
            r_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i),j+2*(N1i+N2i+N3i))       =     -w4i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i),j+N4i+2*(N1i+N2i+N3i))    =    w4i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        
    end
    for j=1:Nwi
        r = [rhowi(j)*cos(thetawi(j)) rhowi(j)*sin(thetawi(j))];
        if j~=Nwi
            r_tmp = [rhowi(j+1)*cos(thetawi(j+1)) rhowi(j+1)*sin(thetawi(j+1))];
        else
            r_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i),j+2*(N1i+N2i+N3i+N4i))       =     -wwi(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i),j+Nwi+2*(N1i+N2i+N3i+N4i))    =    wwi(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        
    end
        b_state(i+N1i+N2i) =  -incfn(norm(p+0.5*p_ed - Rho)) ;

end
for i=1:N4i
    p = [rho4i(i)*cos(theta4i(i)) rho4i(i)*sin(theta4i(i))];
    if i~=N4i
        p_tmp = [rho4i(i+1)*cos(theta4i(i+1)) rho4i(i+1)*sin(theta4i(i+1))];
    else
        p_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1i
        r = [rho1i(j)*cos(theta1i(j)) rho1i(j)*sin(theta1i(j))];
        if j~=N1i
            r_tmp = [rho1i(j+1)*cos(theta1i(j+1)) rho1i(j+1)*sin(theta1i(j+1))];
        else
            r_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i),j)       =     -w1i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i),j+N1i)    =    w1i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        
    end
    for j=1:N2i
        r = [rho2i(j)*cos(theta2i(j)) rho2i(j)*sin(theta2i(j))];
        if j~=N2i
            r_tmp = [rho2i(j+1)*cos(theta2i(j+1)) rho2i(j+1)*sin(theta2i(j+1))];
        else
            r_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i),j+(2*N1i))       =     -w2i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i),j+N2i+(2*N1i))    =    w2i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        
    end
    for j=1:N3i
        r = [rho3i(j)*cos(theta3i(j)) rho3i(j)*sin(theta3i(j))];
        if j~=N3i
            r_tmp = [rho3i(j+1)*cos(theta3i(j+1)) rho3i(j+1)*sin(theta3i(j+1))];
        else
            r_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i),j+2*(N1i+N2i))       =     -w3i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i),j+N3i+2*(N1i+N2i))    =    w3i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        
    end
    for j=1:N4i
        r = [rho4i(j)*cos(theta4i(j)) rho4i(j)*sin(theta4i(j))];
        if j~=N4i
            r_tmp = [rho4i(j+1)*cos(theta4i(j+1)) rho4i(j+1)*sin(theta4i(j+1))];
        else
            r_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A_state(i+(N1i+N2i+N3i),j+2*(N1i+N2i+N3i))       = gdiag4s(k0,i);
            A_state(i+(N1i+N2i+N3i),j+2*(N1i+N2i+N3i)+N4i)    = -1/2.0;
            
        else
            A_state(i+(N1i+N2i+N3i),j+2*(N1i+N2i+N3i))       =  -w4i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A_state(i+(N1i+N2i+N3i),j+2*(N1i+N2i+N3i)+N4i)    = w4i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
            
        end
    end
    for j=1:Nwi
        r = [rhowi(j)*cos(thetawi(j)) rhowi(j)*sin(thetawi(j))];
        if j~=Nwi
            r_tmp = [rhowi(j+1)*cos(thetawi(j+1)) rhowi(j+1)*sin(thetawi(j+1))];
        else
            r_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i),j+2*(N1i+N2i+N3i+N4i))       =     -wwi(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i),j+Nwi+2*(N1i+N2i+N3i+N4i))    =    wwi(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
        
    end
        b_state(i+N1i+N2i+N3i) = - incfn(norm(p+0.5*p_ed - Rho)) ;

end
for i=1:Nwi
    p = [rhowi(i)*cos(thetawi(i)) rhowi(i)*sin(thetawi(i))];
    if i~=Nwi
        p_tmp = [rhowi(i+1)*cos(thetawi(i+1)) rhowi(i+1)*sin(thetawi(i+1))];
    else
        p_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
    end
    p_ed = p_tmp - p;
    for j=1:N1i
        r = [rho1i(j)*cos(theta1i(j)) rho1i(j)*sin(theta1i(j))];
        if j~=N1i
            r_tmp = [rho1i(j+1)*cos(theta1i(j+1)) rho1i(j+1)*sin(theta1i(j+1))];
        else
            r_tmp = [rho1i(1)*cos(theta1i(1)) rho1i(1)*sin(theta1i(1))];
        end
        r_ed = r_tmp - r;
        nhat1 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i+N4i),j)       =     -w1i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i+N4i),j+N1i)    =    w1i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat1),2);
        
    end
    for j=1:N2i
        r = [rho2i(j)*cos(theta2i(j)) rho2i(j)*sin(theta2i(j))];
        if j~=N2i
            r_tmp = [rho2i(j+1)*cos(theta2i(j+1)) rho2i(j+1)*sin(theta2i(j+1))];
        else
            r_tmp = [rho2i(1)*cos(theta2i(1)) rho2i(1)*sin(theta2i(1))];
        end
        r_ed = r_tmp - r;
        nhat2 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i+N4i),j+(2*N1i))       =     -w2i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i+N4i),j+N2i+(2*N1i))    =    w2i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat2),2);
        
    end
    for j=1:N3i
        r = [rho3i(j)*cos(theta3i(j)) rho3i(j)*sin(theta3i(j))];
        if j~=N3i
            r_tmp = [rho3i(j+1)*cos(theta3i(j+1)) rho3i(j+1)*sin(theta3i(j+1))];
        else
            r_tmp = [rho3i(1)*cos(theta3i(1)) rho3i(1)*sin(theta3i(1))];
        end
        r_ed = r_tmp - r;
        nhat3 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i+N4i),j+2*(N1i+N2i))       =     -w3i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i+N4i),j+N3i+2*(N1i+N2i))    =    w3i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat3),2);
        
    end
    for j=1:N4i
        r = [rho4i(j)*cos(theta4i(j)) rho4i(j)*sin(theta4i(j))];
        if j~=N4i
            r_tmp = [rho4i(j+1)*cos(theta4i(j+1)) rho4i(j+1)*sin(theta4i(j+1))];
        else
            r_tmp = [rho4i(1)*cos(theta4i(1)) rho4i(1)*sin(theta4i(1))];
        end
        r_ed = r_tmp - r;
        nhat4 = [r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        
        A_state(i+(N1i+N2i+N3i+N4i),j+2*(N1i+N2i+N3i))       =     -w4i(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
        A_state(i+(N1i+N2i+N3i+N4i),j+N4i+2*(N1i+N2i+N3i))    =    w4i(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhat4),2);
        
    end
    for j=1:Nwi
        r = [rhowi(j)*cos(thetawi(j)) rhowi(j)*sin(thetawi(j))];
        if j~=Nwi
            r_tmp = [rhowi(j+1)*cos(thetawi(j+1)) rhowi(j+1)*sin(thetawi(j+1))];
        else
            r_tmp = [rhowi(1)*cos(thetawi(1)) rhowi(1)*sin(thetawi(1))];
        end
        r_ed = r_tmp - r;
        nhatw = -[r_ed(2)/norm(r_ed) -r_ed(1)/norm(r_ed)];
        if i==j
            A_state(i+(N1i+N2i+N3i+N4i),j+2*(N1i+N2i+N3i+N4i))       = gdiagws(k0,i);
            A_state(i+(N1i+N2i+N3i+N4i),j+2*(N1i+N2i+N3i+N4i)+Nwi)    = -1/2.0;
            
        else
            A_state(i+(N1i+N2i+N3i+N4i),j+2*(N1i+N2i+N3i+N4i))       =  -wwi(j) * glquad(@(t)green(k0,p,0.5,p_ed,r,t,r_ed),2);
            A_state(i+(N1i+N2i+N3i+N4i),j+2*(N1i+N2i+N3i+N4i)+Nwi)    = wwi(j) * glquad(@(t)gradgreen(k0,p,0.5,p_ed,r,t,r_ed,nhatw),2);
            
        end
    end
        b_state(i+N1i+N2i+N3i+N4i) =  -incfn(norm(p+0.5*p_ed - Rho)) ;

end

% disp('State equations populated');

fname = 'state_eqn.mat';
save(fname,'A_state','b_state'); 
