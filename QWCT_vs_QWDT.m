%% Confronto tra il Quantum Walk a tempo discreto e continuo su di una linea infinita
clc;clear;close all;

%il confronto viene fatto tunando gamma in modo da far coincidere
%l'evoluzione a tempo continuo con quella a tempo discreto. La variazione
%di Gamma è fatta per avere un'analoga varianza
tic

%% EVOLUZIONE QWDT

n_step_DT = 100;   %numero di evoluzioni discrete
N = 2*n_step_DT + 3; %serve a garantire l'idealità dell'evoluzione senza subire gli effetti di bordo

N = N_deve_essere_dispari(N);
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

stato_iniziale_DT = kron(sito,coin);
stato_DT = stato_iniziale_DT;


%% Consideriamo ora l'evoluzione
record_evoluzione_DT = zeros (length(stato_DT), n_step_DT + 1); %servirà a salvare il vettore stato ad ogni tempo
record_evoluzione_DT(:,1) = stato_iniziale_DT;


for ii = 1:n_step_DT
    stato_DT = U * stato_DT;
    record_evoluzione_DT(:,ii+1) = stato_DT;
end

%calcolo delle probabilità: ho che il vec stato è del tipo stato = (... jxC;jxT; j+1xC;j+1xT)
%sommo dunque il modulo quadro delle coppie testa croce per avere la
%probabilità di essere al sito j-simo
probabilita_DT = zeros(N,n_step_DT + 1);
varianza_DT = zeros (n_step_DT,1);

for jj = 1:n_step_DT + 1
    for ii = 1:N
        probabilita_DT(ii,jj) = sum(abs(record_evoluzione_DT(2*(ii-1)+1:2*ii,jj)).^2);      
    end
    varianza_DT(jj) = [-(N-1)/2:(N-1)/2].^2*probabilita_DT(:,jj);
end


%% EVOLUZIONE A TEMPO CONTINUO
%numero di siti; numero campionamenti; spaziatura temporale

n_CT_per_DT = 5;%numero di campionamenti in CT rispetto al DT
n_sample_CT = n_CT_per_DT*(n_step_DT -1); %così da avere un andamento più smooth

dt = 1/n_CT_per_DT;

gamma = fit_gamma(varianza_DT);

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
primo_momento_CT = zeros (n_sample_CT,1);
varianza_CT = zeros (n_sample_CT,1);
% Record =dove salviamo l'evoluzione
probabilita_CT = zeros(N,n_sample_CT );
for t = [0:dt:dt*(n_sample_CT-1)]
    ii = ii + 1;
    stato = expm(-1i*L*gamma*t)*stato_iniziale;
    probabilita_CT(1:N,ii) = abs(stato).^2;
    primo_momento_CT(ii) = [1:N]*probabilita_CT(:,ii);
    varianza_CT(ii) = [1:N].^2*probabilita_CT(:,ii) - primo_momento_CT(ii).^2;
end


%% PLOTS
f_50=figure(50);
f_50.Position = [100 200 1500 800];

asse_x = [-(N-1)/2:(N-1)/2];


video = VideoWriter('Video/QWCTvsQWDT.mkv'); % Name it.
video.FrameRate = 10; % How many frames per second.
open(video); 
jj = 1;
for ii = 1:n_sample_CT 

    
    %quando arrivo a tempi interi faccio evolvere il discreto
    if floor ((ii - 1) * dt)== (ii - 1) * dt
        if ii ~= 1
            delete (DT_plot)
        end
        DT_plot = bar(asse_x,probabilita_DT(:,jj), 0.4,'red');  
        hold on
        jj = jj +1;
    end

    if ii ~= 1
        delete (CT_plot)
    end

    CT_plot = plot(asse_x,probabilita_CT(:,ii),'c','LineWidth',2);
    title (['Distribuzione di Probabilità dopo ' , num2str(dt*(ii - 1)) , ' s'])

    xlabel('n^o nodo')
    ylabel('Probabilità')

    xlim ([-(N-1)/2 (N-1)/2])
    if ii < 10
        ylim ([0 1])
    else
        ylim([0 0.15])
    end
    hold on
    legend( 'Tempo Discreto','Tempo Continuo')
    
    pause(.01)
    frame = getframe(gcf); 
    writeVideo(video, frame);

    
  
end
close (video);

toc



function gamma = fit_gamma(varianza_DT)
% la funzione calcola il gamma in modo da fittare l'evoluzione a tempi
% discreti e quella a tempi continui.
%Dalla teoria sappiamo che la varianza è pari a 2*gamma^2*t^2
%Calcoliamo dunque tramite regressione lineare la pendenza di sigma^2 in
%funzione di t^2 (che è uguale a n_step^2) e otteniamo così gamma

%per come è fatta la regressione si veda 
% <https://it.mathworks.com/help/matlab/data_analysis/linear-regression.html> 
t_square = [0:length(varianza_DT) - 1]'.^2;
slope = t_square\varianza_DT;
gamma = sqrt(slope/2);
return

end




function N_out = N_deve_essere_dispari(N)
    if mod(N,2) == 0;
        N_out = N+1;
    else 
        N_out = N;
    end
    return
end
