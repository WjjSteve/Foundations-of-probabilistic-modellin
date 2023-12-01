
M = 1000;
N = 2^15;
T = 30;
dt = T/N;
X0 = 500;   % Initial population
r = .2;   % Intrinsic growth rate
K = 10000;  % Carrying capacity
sigma = .1;  % Intensity of the noise

err = zeros(M, 5);
err2 = zeros(M, 5);

pp = [1 6 7 8 9 10];

for s = 1:M
    dt = T/N;
    dW = sqrt(dt) * randn(1, N);
    W = cumsum(dW);
   

    S1_exa = 0;
    X1_exa = 0;
    S2_exa = 0;
    X2_exa = 0;
    for ind = 1:6
        p = pp(ind);
        R = 2^(p - 1);
        Dt = R * dt;
        L = N/R;
        X1 = X0;
        X2 = X0;
        Y = zeros(1,L);
  
        S1 = 0;
        S2 = 0;
        for j = 1:L
            Winc = sum(dW(R * (j - 1) + 1:R * j));
            X1 = X1 + Dt * r * X1 * (1 - X1/K) + sigma * X1 * Winc;
            X2 = X2  + Dt *r * X2* (1 - X2/K) + sigma * X2 * Winc + 0.5 * X2 *sigma^2 * (Winc.^2-Dt);
            S1 = S1+X1;
            S2 = S2+X2;
        end
        S1 = S1 / L;
        S2 = S2 / L;

        if(ind==1)
            S1_exa = S1;
            S2_exa = S2;
            X1_exa = X1;
            X2_exa = X2;
        else
            err(s, ind-1) = abs(S1 - S1_exa);
            err2(s, ind-1) = abs(S2 - S2_exa);
        end
    end
end
 

size = dt * (2 .^ [5:9]);
er1 = mean(err)
er2 = mean(err2)
%deviation1 = 2*std(err)/sqrt(M);
%deviation2 = 2*std(err2)/sqrt(M);
pl = loglog(size,er1,size,er2)

%errorbar(size,er1, deviation1, 'red',size,er2, deviation2, 'blue');
pl.set('linewidth', 1.3)
set(pl(1), 'marker', 'o')
set(pl(2), 'marker', '^')
dev = std(err2)
ll = legend([pl], {'_D05_';'_D10_'})
xlabel('_xx_');
ylabel('_yy_');

ylim([2e0,1e2])
xlim([2e-2,8e-1])
ll.set('FontSize', 14)
ll.set('Position', [0.6 0.3 0.1 0.2])
ll.set('EdgeColor', 'white')

set(gca,'XTick',[1e-2 1e-1 1e0])
set(gca,'YTick',[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3])
set(gca,'FontSize',14)
