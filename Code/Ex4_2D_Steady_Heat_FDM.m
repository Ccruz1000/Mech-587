% MECH 479 - CFD
% EX4: Two-dimensional steady heat problem by SOR
%------------------------------------------------------------
n=40;
m=40;
iterations=5000;
length=2.0;
h=length/(n-1);
T=zeros(n,m);
bb=1.7;
T(10:n-10,1)=1.0;
for l=1:iterations,
for i=2:n-1, for j=2:m-1
T(i,j)=bb*0.25*(T(i+1,j)+...
T(i,j+1)+T(i-1,j)+T(i,j-1))+(1.0-bb)*T(i,j);
end,end
% find residual
res=0;
for i=2:n-1, 
    for j=2:m-1
        res=res+abs(T(i+1,j)+...
        T(i,j+1)+T(i-1,j)+T(i,j-1)-4*T(i,j))/h^2;
    end
end
l,res/((m-2)*(n-2)) % Print iteration and residual
if (res/((m-2)*(n-2)) < 0.001), 
    break
end
end;
contour(T);