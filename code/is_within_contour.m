function [status] = is_within_contour(p,lambda)
% Returns true if point is within the square bounding the objects
status = true;
% exact
% x1 =  -2.7*lambda; x2 =  -1.3*lambda; y1 =  0.8*lambda; y2 =  2.2*lambda;
% cx  = 2*lambda;   cy = 2*lambda;
% x5 =  -0.2*lambda; x6 =  2.2*lambda; y5 =  -3.075*lambda; y6 =  -1.925*lambda;
% cx1 = -2.5*lambda; cy1= -1*lambda;

% inexact
x1 =  -2.825*lambda; x2 =  -1.175*lambda; y1 =  0.675*lambda; y2 =  2.325*lambda;
x3 =   0.8*lambda; x4 =  3.2*lambda; y3 =  0.8*lambda; y4 =  3.2*lambda;
x5 =  -0.45*lambda; x6 =  2.45*lambda; y5 =  -3.325*lambda; y6 =  -1.675*lambda;
x7 =  -3.95*lambda; x8 =  -1.05*lambda; y7 =  -2.45*lambda; y8 =  0.45*lambda;
% x7 =  -3.8*lambda; x8 =  -1.2*lambda; y7 =  -2.3*lambda; y8 =  0.3*lambda;


if ~(((x1 < p(1))&& (p(1) < x2)) && ((y1 < p(2)) && (p(2) < y2)))
    if ~(((x3 < p(1))&& (p(1) < x4)) && ((y3 < p(2)) && (p(2) < y4)))% Non Exact
%            if ~(sqrt((p(1)-cx)^2 + (p(2)-cy)^2)<=0.85*lambda) % Exact
        if ~(((x5 < p(1))&& (p(1) < x6)) && ((y5 < p(2)) && (p(2) < y6)))
%                          if ~(sqrt((p(1)-cx1)^2 + (p(2)-cy1)^2)<=1.4*lambda)
            if ~(((x7 < p(1))&& (p(1) < x8)) && ((y7< p(2)) && (p(2) < y8)))
                status = false;
                
            end
            
        end
    end
    %    end
end
end

