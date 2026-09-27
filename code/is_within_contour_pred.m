function [status] = is_within_contour_pred(p,lambda)
% Returns true if point is within the square bounding the objects
status = true;
% exact
% x1 =  -2.6*lambda; x2 =  -1.4*lambda; y1 =  0.9*lambda; y2 =  2.1*lambda;
% cx  = 2*lambda;   cy = 2*lambda;
% x5 =  -0.1*lambda; x6 =  2.1*lambda; y5 =  -2.975*lambda; y6 =  -2.025*lambda;
% cx1 = -2.5*lambda; cy1= -1*lambda;


% inexact
x1 =  -2.725*lambda; x2 =  -1.275*lambda; y1 =  0.775*lambda; y2 =  2.225*lambda;
x3 =   0.9*lambda; x4 =  3.1*lambda; y3 =  0.9*lambda; y4 =  3.1*lambda;
x5 =  -0.35*lambda; x6 =  2.35*lambda; y5 =  -3.225*lambda; y6 =  -1.775*lambda;
x7 =  -3.85*lambda; x8 =  -1.15*lambda; y7 =  -2.35*lambda; y8 =  0.35*lambda;
% x7 =  -3.7*lambda; x8 =  -1.3*lambda; y7 =  -2.2*lambda; y8 =  0.2*lambda;

if ~(((x1 < p(1))&& (p(1) < x2)) && ((y1 < p(2)) && (p(2) < y2)))
    if ~(((x3 < p(1))&& (p(1) < x4)) && ((y3 < p(2)) && (p(2) < y4)))% Non Exact
%      if ~(sqrt((p(1)-cx)^2 + (p(2)-cy)^2)<=0.85*lambda) % Exact
        if ~(((x5 < p(1))&& (p(1) < x6)) && ((y5 < p(2)) && (p(2) < y6)))
%             if ~(sqrt((p(1)-cx1)^2 + (p(2)-cy1)^2)<=1.3*lambda) % Exact
             if ~(((x7 < p(1))&& (p(1) < x8)) && ((y7< p(2)) && (p(2) < y8)))
                status = false;
                
            end
            
        end
    end
    %    end
end
end

