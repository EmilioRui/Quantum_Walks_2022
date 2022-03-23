clc;clear;
%% CASO LINEA INFINITA
% Consideriamo DTQW in cui possiamo trascurare effetti di bordo
tic
n_step = 100;   %numero di evoluzioni discrete
N = 2*n_step + 3; %serve a garantire l'idealità dell'evoluzione senza subire gli effetti di bordo

%CONSIDERIAMO L'UNITARIA DI EVOLUZIONE U = S (I x C) 
%S vive nello spazio di Hilbert Hw x Hc, C solo in Hc (coin)

%consideriamo come evoluzione del Coin la matrice di hadamard
C = 1/sqrt(2).*[1 1;1 -1];

%matice con 1 sotto la diagonale, causa j -> j+1
sub_diagonale = zeros(N);
idx = (N+1)* [0:N-2] + 2;
sub_diagonale (idx) = 1;
%matrice con 1 oltre la diagonale causa j -> j-1
over_diagonale = zeros(N);
idx = (N+1)*[1:N-1];
over_diagonale (idx) = 1;

%S causa dunque j->j+1 (jump a dx) se Coin=(1 0) e j->j-1 (sx) se Coin=(0 1)
S = kron(sub_diagonale,[1 0; 0 0]) + kron(over_diagonale,[0 0; 0 1]);

%definiamo dunque U = S (I x C) 
U = S * kron(eye(N),C);

% stato iniziale nel sito 0 (sito centrale della nostra linea),con valore di coin 0
sito = zeros(N,1);
sito((N+1)/2) = 1; %particella inizialmente localizzata
coin = [1 i]'.*(1/sqrt(2)); % stato iniziale di coin 

stato_iniziale = kron(sito,coin);

stato = stato_iniziale;
%% Consideriamo ora l'evoluzione
record_evoluzione = zeros (length(stato), n_step + 1); %servirà a salvare il vettore stato ad ogni tempo
record_evoluzione(:,1) = stato_iniziale;


for ii = 1:n_step
    stato = U * stato;
    record_evoluzione(:,ii+1) = stato;
end

%calcolo delle probabilità: ho che il vec stato è del tipo stato = (... jxC;jxT; j+1xC;j+1xT)
%sommo dunque il modulo quadro delle coppie testa croce per avere la
%probabilità di essere al sito j-simo
probabilita = zeros(N,n_step + 1);
varianza = zeros (n_step,1);
primo_momento = zeros (n_step,1);

for jj = 1:n_step + 1
    for ii = 1:N
        probabilita(ii,jj) = sum(abs(record_evoluzione(2*(ii-1)+1:2*ii,jj)).^2);
        
    end
    primo_momento(jj) = [-(N-1)/2:(N-1)/2]*probabilita(:,jj);
    varianza(jj) = [-(N-1)/2:(N-1)/2].^2*probabilita(:,jj);
end
toc


%% Andiamo in conclusione a plottare i risultati
f_1 = figure(1);
f_1.Position = [100 100 1900 1080];

asse_x = [-(N-1)/2:(N-1)/2]; %e.g N = 5 l'asse va da -2 a 2
% even_mask = mod(asse_x,2)==0;
% uneven_mask = mod(asse_x,2) == 1;

%Dopo un numero pari di jump solo i siti pari hanno probabilità non nulle e
%vicecersa

writerObj = VideoWriter('Video/QWDT_Infinite_line'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 
for ii = 1:n_step + 1    
    bar(asse_x,probabilita(:,ii), 0.3,'red')
    title (['Distribuzione di Probabilità dopo ' , num2str(ii - 1) , ' step'])
    xlim ([-(N-1)/2 (N-1)/2]);
    if ii < 10
        ylim ([0 1])
    else
        ylim([0 0.15])
    end
    pause(.01)
    %prendiamo il video
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);

end
close(writerObj); 



%Controllo che il primo momento venga nullo
figure(2)

plot([0:n_step], primo_momento)
title('primo momento')
ylim ([-0.2 1])
xlabel('n_step')
ylabel('primo momento')

figure(3)
plot([0:n_step].^2, varianza,'r' , 'LineWidth',3)
title('varianza')
xlabel('n_{step}^2')
ylabel ('varianza')



