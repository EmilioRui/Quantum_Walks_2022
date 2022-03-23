clc;clear;close all;
%% Quantum Walk a tempo continuo: Grafo ad anello
%numero di siti; numero campionamenti; spaziatura temporale
N = 100;
n_sample = 200;
dt = 2;

rec_video = false;
%parametri su cui lavorare: gamma ci permette di tunare la velocità di
%propagazione
gamma = 0.2;

N = N_deve_essere_dispari(N);

%definiamo il grado dalla matrice di Adiacenza
A_ring = zeros(N);
% 1 nella sottodiagonale
idx = (N+1)*[0:N-2] + 2;
A_ring (idx) = 1;
% 1 sulla 'sovradiagonale'
idx = (N+1)*[1:N-1];
A_ring(idx) = 1;
%Aggiungiamo le condizioni cicliche
A_ring(1,end) = 1;
A_ring(end,1) = 1;


%definiamo il grafo e ne otteniamo la matrice laplaciana (si potrebbe
%tranquillamente calcolare 'a mano')
Graph_ring = graph(A_ring);
L = laplacian(Graph_ring);


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

f_20 = figure(20);
f_20.Position = [100 100 1500 800];

subplot (1,2,2)
asse_x = [1:N];

if rec_video
    video = VideoWriter('Video/QWCT_ring_graph'); % Name it.
    video.FrameRate = 10; % How many frames per second.
    open(video); 
end
for ii = 1:n_sample + 1
    subplot (1,2,1)
    bar(asse_x,probabilita(:,ii), 0.3,'red')
    title (['Distribuzione di Probabilità dopo ' , num2str((ii - 1)*dt)  , 's'])

    xlabel('n^o nodo')
    ylabel('Probabilità')
    xlim ([1 N])
    if ii < 10
        ylim ([0 1])
    else
        ylim([0 0.15])
    end
    subplot (1,2,2)
    plot(Graph_ring,'EdgeColor','g','NodeColor','c','MarkerSize',100*probabilita(:,ii)+1)
    title(['Grafico ad anello con ', num2str(N), ' siti'])
    
    pause(.03)
    %prendiamo il video
    if rec_video
        frame = getframe(gcf); 
        writeVideo(video, frame);
    end
end
if rec_video 
    close(video);
end
%% Andamento del primo e del secondo momento 

%Controllo che il primo momento venga nullo
figure(22)
plot([0:n_sample], primo_momento)
title('primo momento')
ylim ([N/2-2;N/2+2])
xlabel('n_step')
ylabel('primo momento')

figure(23)
plot([0:n_sample].^2, varianza)
title('varianza')
xlabel('n_step^2')
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


