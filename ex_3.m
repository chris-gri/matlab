load 'Lin_13_sin_pot.mat'


XYZ=XYZ';
R=R';
Q=ElectroStaticBalls(XYZ,R, F);
[f, X, Y, P]=SpherePotential(XYZ, Q, R, r0, a, b, Dx, Dy, Nxy);
figure; 
hold on; 
grid on;
mesh(X,Y,f); 


% [Q, D] = ElectroStaticDipoles(XYZ, R, F);
% Q'
% [f, X, Y, P] = SphereDipPotential(XYZ, Q, D, R, r0, a, b, Dx, Dy, Nxy);
% figure;
% hold on;
% grid on;
% mesh(X,Y,f);