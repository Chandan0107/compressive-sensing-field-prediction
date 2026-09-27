function [P1,P2] = sampling_random(lambda,Nob)
    xl  = -4.55*lambda;
    xh  =  4.55*lambda;
    yl  = -4.55*lambda;
    yh  =  4.55*lambda;
    p1   = randi([xl xh],1,Nob);                                               % generating n points for measurements p1 for x coordinate p2 for y coordinate
    p2   = randi([yl yh],1,Nob);                                               % generating n points for measurements p1 for x coordinate p2 for y coordinate
    x=0;
    for i = 1 : length(p1)
        if ~is_within_contour([p1(i),p2(i)],lambda)
            x = x + 1;
            P1(x) = p1(i);
            P2(x) = p2(i);
        else
            xx = i;    
        end
    end
end











