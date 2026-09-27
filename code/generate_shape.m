function [theta1, rho1, w1, N1] = generate_shape(Zx,Zy,lambda,cx,cy,dis,shape)
% If circle or square, give same Zx and Zy.
    if(shape=='rect')
        x1n1   = -Zx*lambda/2 +cx ; x1n2 = Zx*lambda/2 +cx ;       % X node
        y1n1   = -Zy*lambda/2 +cy ; y1n2 = Zy*lambda/2 +cy ;       % Y node
        Ny     = ceil(dis*((y1n2 - y1n1)/lambda));                   % Number of segments   
        Nx     = ceil(dis*((x1n2 - x1n1)/lambda));
        %Ndis1  = N1/4;
        disty  = (y1n2 - y1n1)/Ny;
        distx  = (x1n2 - x1n1)/Nx;
        xcor1  = (x1n1:distx:x1n2-distx);
        ycor1  = y1n1*ones(1,Nx);
        xcor2  = x1n2*ones(1,Ny);
        ycor2  = (y1n1:disty:y1n2-disty);
        xcor3  = (x1n2:-distx:x1n1+distx);
        ycor3  = y1n2*ones(1,Nx);
        xcor4  = x1n1*ones(1,Ny);
        ycor4  = (y1n2:-disty:y1n1+disty);
        x1     = [xcor1 xcor2 xcor3 xcor4];
        y1     = [ycor1 ycor2 ycor3 ycor4];
        n1     = max(size(x1));
        [theta1,rho1] = cart2pol(x1,y1);
        N1=2*(Nx+Ny);
    elseif(shape=='circ')
        r = Zx*lambda; % radius of the circle
        N1 = round(2*pi*r*dis/lambda);
        th = 0:2*pi/N1:2*pi-2*pi/N1;
        x1 = cx + r*cos(th); y1 = cy + r*sin(th); 
        [theta1,rho1] = cart2pol(x1,y1);
    end
    % Dist b/w points
    w1=zeros(N1,1);
    for i=1:N1-1
        w1(i)=sqrt((rho1(i+1)*cos(theta1(i+1))-rho1(i)*cos(theta1(i)))^2+(rho1(i+1)*sin(theta1(i+1))-rho1(i)*sin(theta1(i)))^2);
    end
    w1(N1)=sqrt((rho1(end)*cos(theta1(end))-rho1(1)*cos(theta1(1)))^2+(rho1(end)*sin(theta1(end))-rho1(1)*sin(theta1(1)))^2);
end

