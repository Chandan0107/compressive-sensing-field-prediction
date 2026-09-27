 function [ y ] = gradgreen(k,p,s,p_ed,r,t,r_ed,n)
%Grad(greenfn).hat(n) y = jk/4 H_1^(2)(k|r-r'|)hat(r-r') dot hat(n)
ps = p+s*p_ed;
rt = r+t*r_ed;
rp = rt - ps; 
rho = norm(rp);
y = 1j*k/4*besselh(1,2,k*rho)*dot(rp,n)/rho;
end
