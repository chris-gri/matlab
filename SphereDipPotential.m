function [F,X,Y,P] = SphereDipPotential(XYZ,Q,D, R,r0,a,b,Dx,Dy,Nxy)
P = zeros(3,2);
c=b-(b'*a)/(a'*a)*a;
if norm(c)<10^(-10)
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
      X(i, j)=Dx(1) + (j - 1)/(Nxy(2) - 1)*(Dx(2) - Dx(1));
      Y(i, j)=Dy(1) + (i - 1)/(Nxy(1) - 1)*(Dy(2) - Dy(1));
      r1 = r0 + P*[X(i, j);Y(i, j)];
      for k = 1:1:n
          r = r1 - XYZ(:, k);
          d = norm(r);
          if d > R(k)
             F(i, j) = F(i, j) + (D(k, :)*r)/(d^3) + Q(k)/d;
          else
             F(i, j) = F(i, j) + (D(k, :)*r)/(R(k)^3) + Q(k)/R(k);
          end
      end
   end
end
end
