function [Q, D] = ElectroStaticDipoles(XYZ,R,F)
    R'
    N = length(R)
    A = zeros(N)
   
    
    for i = [1:N]
        for j = [1:N]
            if i==j
                A(i, i) = 1/R(i);
            else 
                A(i, j) = ((XYZ(i,1) - XYZ(j,1))^2 + (XYZ(i,2) - XYZ(j,2))^2 + (XYZ(i,3) - XYZ(j,3))^2)^(-1/2);
                if (R(i)+R(j))*A(i, j) >= 1
                    error('too close')
                end
            end
        end
    end
    
    B = zeros(4*N);
    for i = [1:N]
        for j = [1:N]
            B(i,j) = A(i,j)
            B(i,j+N) =   (XYZ(i,1) - XYZ(j,1))*A(i,j)^3
            B(i,j+2*N) = (XYZ(i,2) - XYZ(j,2))*A(i,j)^3
            B(i,j+3*N) = (XYZ(i,3) - XYZ(j,3))*A(i,j)^3
            
            B(i+N,j) = (XYZ(i,1) - XYZ(j,1))*A(i,j)^3
            B(i+N,j+N) =   3*(XYZ(i,1) - XYZ(j,1))*(XYZ(i,1) - XYZ(j,1))*A(i,j)^5 - A(i,j)^3
            B(i+N,j+2*N) = 3*(XYZ(i,1) - XYZ(j,1))*(XYZ(i,2) - XYZ(j,2))*A(i,j)^5
            B(i+N,j+3*N) = 3*(XYZ(i,1) - XYZ(j,1))*(XYZ(i,3) - XYZ(j,3))*A(i,j)^5
            
            B(i+2*N,j) = (XYZ(i,2) - XYZ(j,2))*A(i,j)^3
            B(i+2*N,j+N) =   3*(XYZ(i,2) - XYZ(j,2))*(XYZ(i,1) - XYZ(j,1))*A(i,j)^5
            B(i+2*N,j+2*N) = 3*(XYZ(i,2) - XYZ(j,2))*(XYZ(i,2) - XYZ(j,2))*A(i,j)^5 - A(i,j)^3
            B(i+2*N,j+3*N) = 3*(XYZ(i,2) - XYZ(j,2))*(XYZ(i,3) - XYZ(j,3))*A(i,j)^5
            
            B(i+3*N,j) = (XYZ(i,3) - XYZ(j,3))*A(i,j)^3
            B(i+3*N,j+N) =   3*(XYZ(i,3) - XYZ(j,3))*(XYZ(i,1) - XYZ(j,1))*A(i,j)^5
            B(i+3*N,j+2*N) = 3*(XYZ(i,3) - XYZ(j,3))*(XYZ(i,2) - XYZ(j,2))*A(i,j)^5
            B(i+3*N,j+3*N) = 3*(XYZ(i,3) - XYZ(j,3))*(XYZ(i,3) - XYZ(j,3))*A(i,j)^5 - A(i,j)^3
        end
    end
    
    E = zeros(3*N,1)
    FE = [F; E]
    QP = B\FE
    Q = QP(1:N,1)'
    D = [QP(N+1:2*N,1), QP(2*N+1:3*N,1), QP(3*N+1:4*N,1)]
end
