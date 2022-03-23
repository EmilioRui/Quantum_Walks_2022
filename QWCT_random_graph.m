clc;clear;
%% Quantum Walk a tempo continuo: Grafo random
%numero di siti; numero campionamenti; spaziatura temporale
N = 100;
n_sample = 100;
dt = 2;
p_link = 5/N; %probabilità di avere il link tra 2 nodi
%parametri su cui lavorare: gamma ci permette di tunare la velocità di
%propagazione
gamma = 0.2;

N = N_deve_essere_dispari(N);



%definiamo il grafo e ne otteniamo la matrice laplaciana (si potrebbe
%tranquillamente calcolare 'a mano')
G = randgraph(N,p_link);
L = laplacian(G);


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

f_40 = figure(40)
f_40.Position = [100 100 1500 800];

subplot (1,2,2)
asse_x = [1:N];
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
    plot(G,'EdgeColor','g','NodeColor','r','MarkerSize',50*probabilita(:,ii)+1)
    title(['Grafico ad anello con ', num2str(N), ' siti'])

    pause(.01)
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






%% Generatore casuale di grafi
function G = randgraph(N,p_1)
%data la dimensione N la funzione genera una matrice di adiacenza casuale e
%da questa un grafo. L'output G è un grafo. p_0/è la frequenza con cui vogliamo le connessioni

    %caso base: probabilità uguale di avere/non avere il link
    if nargin == 1
        p_1 = .5;
    end
    A = rand(N); %generiamo una matrice casuale di numeri compresi tra 0 e 1
    A(A>(1-p_1)) = 1;
    A(A ~= 1)= 0;
    
    %Per ottenere una matrice simmetrica con diagonale nulla prendiamo la parte
    %triangolare superiore, la riportiamo trasponendola nella parte inferiore
    %ed eliminiamo la diagonale
    A = triu(A) + triu(A,1)'- diag(diag(A));
    G = graph(A);
    return
    end