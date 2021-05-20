function [F,X,Y,P] = SpherePotential(XYZ, Q, R, r0, a, b, Dx, Dy, Nxy)
 
P=zeros(2);
c=b-(b'*a)/(a'*a)*a;
if norm(c)<10^(-5)
   error('no way'); 
end
P(:, 1)=a;
P(:, 2)=c;
F=zeros(Nxy(1), Nxy(2));
Y=zeros(Nxy(1), Nxy(2));
X=zeros(Nxy(1), Nxy(2));

n=size(Q);
n=n(1);
for i=1:1:Nxy(1)
   for j=1:1:Nxy(2)
      X(i, j)=Dx(1)+(j-1)/(Nxy(2)-1)*(Dx(2)-Dx(1));
      Y(i, j)=Dy(1)+(i-1)/(Nxy(1)-1)*(Dy(2)-Dy(1));
      r=r0+P*[X(i, j);Y(i, j)];
      for k = 1:1:n
          d = norm(XYZ(:, k)-r);
          if d > R(k)
             F(i, j)=F(i, j) + Q(k)/d;
          else
             F(i, j)=F(i, j) + Q(k)/R(k);
          end
      end
   end
end
end