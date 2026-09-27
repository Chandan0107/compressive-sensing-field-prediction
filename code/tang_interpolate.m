clear all
close all
clc
% load('Four objects/Closed_Test_4obj.mat','xcoe1','lambda','N1','N2','N3','N4','Nw')
load('true_tangfields_all_surfaces_40.mat','tang_fields');
lambda = 20;
xcoe1 = tang_fields;
dis = 40;

% inexact
Zx=1.25; Zy=1.25; cx=-2*lambda; cy=2.5*lambda;  shape='rect';
[theta1, rho1, w1, N1] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape='rect';
[theta2, rho2, w2, N2] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=1.25; cx=2*lambda; cy=-2.5*lambda; shape='rect'; 
[theta3, rho3, w3, N3] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-2*lambda; 
[theta4, rho4, w4, N4] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  
[thetaw, rhow, ww, Nw] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape);

dis1 = 5;
% rect
x_gp1 = xcoe1(1:N1);
x_p1  = xcoe1(N1+1:2*N1);

% circ
x_gp2 = xcoe1(2*N1+1:2*N1+N2);
x_p2  = xcoe1(2*N1+N2+1:2*N1+2*N2);

% rect 
x_gp3  = xcoe1(2*N1+2*N2+1 : 2*(N1+N2)+N3);
x_p3   = xcoe1(2*N1+2*N2+N3+ 1 : 2*(N1+N2+N3));

x_gp4 = xcoe1(2*(N1+N2+N3)+1:2*(N1+N2+N3)+N4);
x_p4  = xcoe1(2*(N1+N2+N3)+N4+1:2*(N1+N2+N3+N4));

% rect
x_gpw = xcoe1(2*(N1+N2+N3+N4)+1:2*(N1+N2+N3+N4)+Nw);
x_pw  = xcoe1(2*(N1+N2+N3+N4)+ Nw +1:end);




%% Object 1
Zx = 1.25; Zy=1.25; cx=-2*lambda; cy=2.5*lambda; shape='rect'; 
[xh1 yh2 xh3 yh4 Nhx Nhy Nh1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis);

Zx = 1.25; Zy=1.25; cx=-2*lambda; cy=2.5*lambda; shape='rect'; 
[xl1 yl2 xl3 yl4 Nlx Nly Nl1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis1);


xl_gp1_1 = interp1(xh1,x_gp1(1:length(xh1)),xl1,'spline');
xl_gp1_2 = interp1(yh2,x_gp1(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_gp1_3 = interp1(xh3,x_gp1(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_gp1_4 = interp1(yh4,x_gp1(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_gp1 = [xl_gp1_1, xl_gp1_2, xl_gp1_3, xl_gp1_4].';

xl_p1_1 = interp1(xh1,x_p1(1:length(xh1)),xl1,'spline');
xl_p1_2 = interp1(yh2,x_p1(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_p1_3 = interp1(xh3,x_p1(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_p1_4 = interp1(yh4,x_p1(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_p1 = [xl_p1_1, xl_p1_2, xl_p1_3, xl_p1_4].';

%% Object 2
Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape='rect';
[xh1 yh2 xh3 yh4 Nhx Nhy Nh1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis);

Zx=2; Zy=2; cx=2*lambda; cy=2*lambda; shape='rect';
[xl1 yl2 xl3 yl4 Nlx Nly Nl1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis1);


xl_gp2_1 = interp1(xh1,x_gp2(1:length(xh1)),xl1,'spline');
xl_gp2_2 = interp1(yh2,x_gp2(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_gp2_3 = interp1(xh3,x_gp2(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_gp2_4 = interp1(yh4,x_gp2(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_gp2 = [xl_gp2_1, xl_gp2_2, xl_gp2_3, xl_gp2_4].';

xl_p2_1 = interp1(xh1,x_p2(1:length(xh1)),xl1,'spline');
xl_p2_2 = interp1(yh2,x_p2(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_p2_3 = interp1(xh3,x_p2(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_p2_4 = interp1(yh4,x_p2(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_p2 = [xl_p2_1, xl_p2_2, xl_p2_3, xl_p2_4].';

%% Object 3
Zx=2.5; Zy=1.25; cx=2*lambda; cy=-2.5*lambda; shape='rect';
[xh1 yh2 xh3 yh4 Nhx Nhy Nh1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis);

Zx=2.5; Zy=1.25; cx=2*lambda; cy=-2.5*lambda; shape='rect';
[xl1 yl2 xl3 yl4 Nlx Nly Nl1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis1);


xl_gp3_1 = interp1(xh1,x_gp3(1:length(xh1)),xl1,'spline');
xl_gp3_2 = interp1(yh2,x_gp3(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_gp3_3 = interp1(xh3,x_gp3(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_gp3_4 = interp1(yh4,x_gp3(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_gp3 = [xl_gp3_1, xl_gp3_2, xl_gp3_3, xl_gp3_4].';

xl_p3_1 = interp1(xh1,x_p3(1:length(xh1)),xl1,'spline');
xl_p3_2 = interp1(yh2,x_p3(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_p3_3 = interp1(xh3,x_p3(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_p3_4 = interp1(yh4,x_p3(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_p3 = [xl_p3_1, xl_p3_2, xl_p3_3, xl_p3_4].';



%% Object 4
Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-2*lambda; shape='rect';
[xh1 yh2 xh3 yh4 Nhx Nhy Nh1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis);

Zx=2.5; Zy=2.5; cx=-2.5*lambda; cy=-2*lambda; shape='rect';
[xl1 yl2 xl3 yl4 Nlx Nly Nl1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis1);


xl_gp4_1 = interp1(xh1,x_gp4(1:length(xh1)),xl1,'spline');
xl_gp4_2 = interp1(yh2,x_gp4(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_gp4_3 = interp1(xh3,x_gp4(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_gp4_4 = interp1(yh4,x_gp4(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_gp4 = [xl_gp4_1, xl_gp4_2, xl_gp4_3, xl_gp4_4].';

xl_p4_1 = interp1(xh1,x_p4(1:length(xh1)),xl1,'spline');
xl_p4_2 = interp1(yh2,x_p4(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_p4_3 = interp1(xh3,x_p4(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_p4_4 = interp1(yh4,x_p4(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_p4 = [xl_p4_1, xl_p4_2, xl_p4_3, xl_p4_4].';




%% Wall
Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  
[xh1 yh2 xh3 yh4 Nhx Nhy Nh1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis);

Zx = 9.5; Zy = 9.5; cx = 0; cy = 0;  
[xl1 yl2 xl3 yl4 Nlx Nly Nl1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis1);


xl_gpw_1 = interp1(xh1,x_gpw(1:length(xh1)),xl1,'spline');
xl_gpw_2 = interp1(yh2,x_gpw(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_gpw_3 = interp1(xh3,x_gpw(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_gpw_4 = interp1(yh4,x_gpw(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_gpw = [xl_gpw_1, xl_gpw_2, xl_gpw_3, xl_gpw_4].';

xl_pw_1 = interp1(xh1,x_pw(1:length(xh1)),xl1,'spline');
xl_pw_2 = interp1(yh2,x_pw(length(xh1)+1:length(xh1)+length(yh2)),yl2,'spline');
xl_pw_3 = interp1(xh3,x_pw(length(xh1)+length(yh2)+1:length(xh1)+length(yh2)+length(xh3)),xl3,'spline');
xl_pw_4 = interp1(yh4,x_pw(length(xh1)+length(yh2)+length(xh3)+1:end),yl4,'spline');

xl_pw = [xl_pw_1, xl_pw_2, xl_pw_3, xl_pw_4].';

tang_fields_l = [xl_gp1; xl_p1; xl_gp2; xl_p2; xl_gp3; xl_p3; xl_gp4; xl_p4; xl_gpw; xl_pw];

figure; 
subplot(211); plot(abs(xcoe1)); title('dis = 40')
subplot(212); plot(abs(tang_fields_l)); title('lower discretization');

fname = strcat('true_tang_fields_dis_',num2str(dis1),'.mat');
save(fname,'tang_fields_l','dis1');

function [xcor1 ycor2 xcor3 ycor4 Nx Ny N1] = generate_coordinates(Zx,Zy,lambda,cx,cy,dis)
% If circle or square, give same Zx and Zy.
        x1n1   = -Zx*lambda/2 +cx ; x1n2 = Zx*lambda/2 +cx ;       % X node
        y1n1   = -Zy*lambda/2 +cy ; y1n2 = Zy*lambda/2 +cy ;       % Y node
        Ny     = ceil(dis*((y1n2 - y1n1)/lambda));                   % Number of segments   
        Nx     = ceil(dis*((x1n2 - x1n1)/lambda));
        %Ndis1  = N1/4;
        disty  = (y1n2 - y1n1)/Ny;
        distx  = (x1n2 - x1n1)/Nx;
        xcor1a  = (x1n1:distx:x1n2-distx);
        xcor1 = [];
        for i = 1:length(xcor1a)-1
            xcor1 = [xcor1, (xcor1a(i) + xcor1a(i+1))/2];
        end
        xcor1 = [xcor1, (xcor1a(end) + x1n2)/2];
        ycor1  = y1n1*ones(1,Nx);
        xcor2  = x1n2*ones(1,Ny);
        ycor2a  = (y1n1:disty:y1n2-disty);
        ycor2 = [];
        for i = 1:length(ycor2a)-1
            ycor2 = [ycor2, (ycor2a(i) + ycor2a(i+1))/2];
        end
        ycor2 = [ycor2, (ycor2a(end) + y1n2)/2];
        xcor3a  = (x1n2:-distx:x1n1+distx);
        xcor3 = [];
        for i = 1:length(xcor3a)-1
            xcor3 = [xcor3, (xcor3a(i) + xcor3a(i+1))/2];
        end
        xcor3 = [xcor3, (xcor3a(end) + x1n1)/2];
        ycor3  = y1n2*ones(1,Nx);
        xcor4  = x1n1*ones(1,Ny);
        ycor4a  = (y1n2:-disty:y1n1+disty);
        ycor4 = [];
        for i = 1:length(ycor4a)-1
            ycor4 = [ycor4, (ycor4a(i) + ycor4a(i+1))/2];
        end
        ycor4 = [ycor4, (ycor4a(end)+y1n1)/2];
        x1     = [xcor1 xcor2 xcor3 xcor4];
        y1     = [ycor1 ycor2 ycor3 ycor4];
        n1     = max(size(x1));
%         [theta1,rho1] = cart2pol(x1,y1);
        N1=2*(Nx+Ny);
end

