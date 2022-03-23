clc;clear;
%% Quantum Walk a tempo continuo: Grafo a Linea
%numero di siti; numero campionamenti; spaziatura temporale
N = 100;
n_sample = 200;
dt = 2;

%decidiamo se riprendere i plot o no
rec_video = false;
%parametri su cui lavorare: gamma ci permette di tunare la velocità di
%propagazione
gamma = 0.2;

N = N_deve_essere_dispari(N);

%definiamo il grado dalla matrice di Adiacenza
A_line = zeros(N);
% 1 nella sottodiagonale
idx = (N+1)*[0:N-2] + 2;
A_line (idx) = 1;
% 1 sulla 'sovradiagonale'
idx = (N+1)*[1:N-1];
A_line(idx) = 1;
%così definita causa già rimbalzi ai bordi


%definiamo il grafo e ne otteniamo la matrice laplaciana (si potrebbe
%tranquillamente calcolare 'a mano')
Graph_line = graph(A_line);
L = laplacian(Graph_line);


%definizione stato iniziale: stato localizzato al centro della linea
stato_iniziale = zeros(N,1);
stato_iniziale((N+1)/2) = 1;

%% Evoluzione
ii = 0;
primo_momento = zeros (n_sample,1);
varianza = zeros (n_sample,1);
% Record =dove salviamo l'evoluzione
probabilita = zeros(N,n_sample +1);
for t = [0:dt:dt*n_sample]
    ii = ii + 1;
    stato = expm(-1i*L*gamma*t)*stato_iniziale;
    probabilita(1:N,ii) = abs(stato).^2;
    primo_momento(ii) = [1:N]*probabilita(:,ii);
    varianza(ii) = [1:N].^2*probabilita(:,ii) - primo_momento(ii).^2;
end

%% Plots



f_10 = figure(10);
f_10.Position = [100 100 1500 800];
subplot (1,2,2)
asse_x = [1:N];
if rec_video
    video = VideoWriter('Video/QWCT_Graph_line'); % Name it.
    video.FrameRate = 10; % How many frames per second.
    open(video); 
end
for ii = 1:n_sample + 1
    subplot (1,2,1)
    bar(asse_x,probabilita(:,ii), 0.3,'red')
    title (['Distribuzione di Probabilità dopo ' , num2str((ii - 1)*dt)  , 's'])
    xlim ([1 N])
    xlabel('n^o nodo')
    ylabel('Probabilità')
    if ii < 10
        ylim ([0 1])
    else
        ylim([0 0.15])
    end
    subplot (1,2,2)
    plot(Graph_line,'EdgeColor','g','NodeColor','r','MarkerSize',100*probabilita(:,ii)+1)
    pause(.01)

    %prendiamo il video
    if rec_video
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(video, frame);
    end
end
if rec_video
    close(video); 
end

%% Andamento del primo e del secondo momento 

%Controllo che il primo momento venga nullo
figure(12)
plot([0:n_sample], primo_momento , 'r')
title('primo momento')
ylim ([N/2-2;N/2+2])
xlabel('n_{step}')
ylabel('primo momento')

figure(13)
plot([0:n_sample].^2, varianza , 'r' , 'LineWidth',3)
title('varianza QWCT grafo a linea')
xlabel('n_{step}^2')
ylabel ('varianza')





%% Per come è stato scritto il programma è comodo lavorare con uno stato iniziale localizzato nel centro 
function N_out = N_deve_essere_dispari(N)
    if mod(N,2) == 0;
        N_out = N+1;
    else 
        N_out = N;
    end
    return
end

