randn('state', 100);



M = 1000;
mu = 5;
sigma = 2; 
N = 2^12;
T = 1;
dt = T/N;
X0 = 10;  


err = zeros(M, 5);

for s = 1:M
    dt = T/N;
    dW = sqrt(dt) * randn(1, N);

    W = cumsum(dW);
    X_exa = X0 * exp((mu - 0.5 * sigma^2) + sigma * W(end));
   
    for p = 1:5  
        R = 2^(p - 1);
        Dt = R * dt;
        L = N/R;
        X = X0;
        Y = zeros(1,L);
        S = 0;
        TE = 0;
        for j = 1:L
            Winc = sum(dW(R * (j - 1) + 1:R * j));
            X = X + Dt * mu * X + sigma * X * Winc;
            %X = X + Dt * mu * X + sigma * X * Winc + 0.5 * X * sigma^2 * (Winc.^2-Dt);
        end
    
      err(s, p) = abs(S - S_exa);
    end
end
 

size = dt * (2 .^ [0:4]);
er = mean(err)
deviation = 2*std(err)/sqrt(M);
pl = loglog(size,size.^.5,size,er)
hold on 
errorbar(size, er, deviation, 'red');
pl.set('linewidth', 1.3)
set(pl(1), 'marker', 'o')
set(pl(2), 'marker')
ll = legend([pl], {'_D05_';'_D10_'})
xlabel('_xx_');
ylabel('_yy_');

ylim([1e-5,1e-2])
xlim([2e-5,8e-4])
ll.set('FontSize', 14)
ll.set('Position', [0.7 0.2 0.1 0.2])
ll.set('EdgeColor', 'white')
