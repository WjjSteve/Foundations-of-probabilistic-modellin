N = 100;        
J = 1;          
      

spins = ones(N,N);
%spins = randi([0, 1], N, N) * 2 - 1;


n = 1000;
S = N * N;

% Initialize magnetization array to store results
magnetization = zeros(1, n);
tau = [];
B = [];
E = [];

% imagesc(spins);
% axis off;
% colormap([1 1 1; 0 0 0]); 
% caxis([-1, 1]);
% title(['T = 2.4, steps = ',num2str(0)]);
% exportgraphics(gca,"parabola.gif","Append",true)

for T = 2.4
    beta = 1 / (J * T); 

    for step = 1:n
        for swip_step = 1:S
            i = randi(N);
            j = randi(N);
            
            delta_E = 2 * J * spins(i, j) * (spins(mod(i+1-1, N)+1, j) + ...
                                             spins(mod(i-1-1, N)+1, j) + ...
                                             spins(i, mod(j+1-1, N)+1) + ...
                                             spins(i, mod(j-1-1, N)+1));
    
          if ((delta_E <= 0) || (rand() <= exp(-beta * delta_E)))
              if(step == 1)
                % imagesc(spins);
                % axis off;
                % colormap([1 1 1; 0 0 0]); 
                % caxis([-1, 1]);
                % title(['T = 2.4, steps = ',num2str(0)]);
                % exportgraphics(gca,"parabola.gif","Append",true)
                spins(i, j) = -spins(i, j);
              end
          end

        end

          % imagesc(spins);
          % axis off
          % colormap([1 1 1; 0 0 0]); 
          % title(['T = 2.4, steps = ',num2str(step*S-S+swip_step)]);
          % exportgraphics(gca,"parabola.gif","Append",true)
   

        magnetization(step) = abs(sum(sum(spins))) / S;
    end

    stable = magnetization(8001:end); 
    [y,lags] = autocorr(stable,'numlags',1999);
    y = y - exp(-1);
    Q = find(y<=0);
    Q = Q(1);
    tau = [tau,Q];
    
    stable_new = magnetization(8001:2*Q:end);
     
    bb = mean(stable_new);
    ee = 2*std(stable_new)/sqrt(size(stable_new,2));
     
    B = [B,bb];
    E = [E,ee];
end
