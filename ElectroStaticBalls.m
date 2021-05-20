function [Q] = ElectroStaticBalls(XYZ, R, F)
s = size(R); 
s = s(2); % second component of the matrix size ps nomber of balls 
C = zeros(s);

for i=1:1:s
   for j=1:1:s
      if i==j
         C(i, j)=1/R(i);
      else
          d=((XYZ(1, i)-XYZ(1, j)).^2+(XYZ(2, i)-XYZ(2, j)).^2+(XYZ(3, i)-XYZ(3, j)).^2).^0.5; % distance between centres of balls 
          if d > R(i) + R(j)
            C(i, j) = 1/d;
          else
            error('too close');
          end
      end
   end
end

Q=C\F;

end