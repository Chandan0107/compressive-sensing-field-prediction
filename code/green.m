function [ y ] = green(k,p,s,p_ed,r,t,r_ed)
%Green's fn y = -j/4 H_0^(2)(k|(p+s*p_ed)-(r+t*r_ed))
ps = p+s*p_ed;
rt = r+t*r_ed;
rp = rt - ps;
y = -1j/4*besselh(0,2,k*norm(rp));
end


