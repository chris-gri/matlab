function [matrix] = KronigPenney(k_input, m_input, a_input, b_input, U0_input, Emax_input)

global k 
global U0
global h
global m
global a
global b
global Emax

h=1.054571817*10.^(-27)/((1.602176634*10.^(-12)*1*10.^(-14)).^0.5);

a = a_input;
b = b_input;
U0 = U0_input;
Emax = Emax_input; 
m = m_input;
f=@F;

final_matrix =[];
matrix_length = 0;

for ii = 1:length(k_input)
    k = k_input(ii);
    x =[];
    for jj = U0+0.0001:0.1:Emax
        x1 = fzero(f, jj);
        y = true;
        for kk = 1:length(x)
            if abs(x(kk) - x1) < 0.0001
                y = false;
            end
        end
        if (y) && (x1 < Emax)
            x = [x ; x1];
        end
    end
    x = sort(x);
    if ii == 1
        final_matrix = x;
        matrix_length = length(x);
    else 
        if matrix_length < length(x)
            final_matrix = [final_matrix; NaN.*zeros(length(x)-matrix_length, ii-1)]; 
            matrix_length = length(x);
        else
            x = [x; NaN.*zeros(matrix_length - length(x), 1)];
        end
        
        final_matrix = [final_matrix, x];
    end
end

matrix = final_matrix;
end

function cos1 = F(E)
global k 
global U0
global h
global m
global a
global b
mu = (2*m.*E/h.^2).^0.5;
lambda = (2*m.*(E-U0)/h.^2).^0.5;

cos1 = cos(mu.*a).*cos(lambda.*b)-(mu.^2+lambda.^2)./(2*mu.*lambda).*sin(mu.*a).*sin(lambda.*b)-cos(k.*(a+b));
end