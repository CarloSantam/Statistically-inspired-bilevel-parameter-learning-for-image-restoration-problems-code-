function fourier_matrix = generate_fourier_matrix(n)
    % Crea una matrice vuota con dimensioni nxn
    fourier_matrix = zeros(n, n);
    
    % Calcola gli elementi della matrice di Fourier
    for i = 1:n
        for j = 1:n
            exponent = 2i * pi * (i-1) * (j-1) / n;
            fourier_matrix(i, j) = exp(exponent) / sqrt(n);
        end
    end
end



%Questo genererà una matrice di Fourier di dimensione 4x4 e la visualizzerà nella finestra di output di MATLAB. Tieni presente che la matrice di Fourier è complessa, quindi gli elementi potrebbero avere parti reali e immaginarie.
