%% ANN
clear
close all
clc

% CARICAMENTO DEI DATI

load -ascii training_set_contest.txt
load -ascii test_set_contest.txt

T = 365; f = 10; 
anni=size(training_set_contest,1)/T;

data = training_set_contest(1:anni*T,2:14); % da 2, tolgo lo storage

% DESTAGIONALIZZAZIONE (su tutto il dataset)

for i = 1:size(data,2)
    [ x(:,i),m(:,i),sd(:,i) ] = destagionalizzo(data(:,i), T, f);
end
% x = dati destagionalizzati
% m = medie mobili ripetute per tutta la serie temporale
% sd = deviazioni standard ripetute per tutta la serie temporale



%Dividiamo il dataset in calibrazione e validazione

ac=ceil(0.7*anni); %ac=anni di calibrazione

dc=data(1:ac*T,:); %dc = dati calibrazione
dv=data(ac*T+1:end,:); %dv = dati validazione

xc=x(1:ac*T,:); %xc = dati calibrazione destagionalizzati
xv=x(ac*T+1:end,:); %xv = dati validazione destagionalizzati

mc = m(1:ac*T,:);
mv = m(ac*T+1:end,:);

sdc = sd(1:ac*T,:);
sdv = sd(ac*T+1:end,:);

%% Partizione del dataset per la cross validazione

%K-fold cross validation

K=anni; %K=anni

ac=ceil(0.7*anni); %ac=anni di calibrazione
av=anni-ac; %av=anni di validazione
num_variabili = size(data,2);

%Creo matrici 3d a K strati, dove ogni strato k sarà un dataset diverso.
%A ogni strato k, gli anni di calibrazione/validazione shiftano di uno
%Ad es. k=1 --> anni calibraz = 1:5, anni validaz =6,7. k=2 --> anni
%calibraz = 2:6, anni validaz = 1,7; etc...

DC = nan(ac*T, num_variabili, K);
DV = nan(av*T, num_variabili, K);

XC = nan(ac*T, num_variabili, K);
XV = nan(av*T, num_variabili, K);

MC = nan(ac*T, num_variabili, K);
MV = nan(av*T, num_variabili, K);

SDC = nan(ac*T, num_variabili, K);
SDV = nan(av*T, num_variabili, K);

%Riempimento con un ciclo for

for k = 1:K
    
    if ac+k-1<=anni
       
        DC(:,:,k)=data( (k-1)*T+1:(ac+k-1)*T , : );
        XC(:,:,k)=x( (k-1)*T+1:(ac+k-1)*T , : );
        MC(:,:,k)=m( (k-1)*T+1:(ac+k-1)*T , : );
        SDC(:,:,k)=sd( (k-1)*T+1:(ac+k-1)*T , : );

    else
        DC(:,:,k)=[data(1:(k-av-1)*T, : ) ; data((k-1)*T+1:anni*T , : )];
        XC(:,:,k)=[x(1:(k-av-1)*T, : ) ; x((k-1)*T+1:anni*T , : )];
        MC(:,:,k)=[m(1:(k-av-1)*T, : ) ; m((k-1)*T+1:anni*T , : )];
        SDC(:,:,k)=[sd(1:(k-av-1)*T, : ) ; sd((k-1)*T+1:anni*T , : )];
    end
        
    if k==1
    
        DV(:,:,k)=data( ac*T+1:anni*T , : );
        XV(:,:,k)=x( ac*T+1:anni*T , : );
        MV(:,:,k)=m( ac*T+1:anni*T , : );
        SDV(:,:,k)=sd( ac*T+1:anni*T , : );
    
    end
    
    if k>av
        
        DV(:,:,k)=data( (k-av-1)*T+1:(k-1)*T , : );
        XV(:,:,k)=x( (k-av-1)*T+1:(k-1)*T , : );
        MV(:,:,k)=m( (k-av-1)*T+1:(k-1)*T , : );
        SDV(:,:,k)=sd( (k-av-1)*T+1:(k-1)*T , : );
        
    end
    
    if k~=1 && k<=av
        
        DV(:,:,k)=[data( 1:(k-1)*T , : ) ; data( (ac+k-1)*T+1:anni*T , : )];
        XV(:,:,k)=[x( 1:(k-1)*T , : ) ; x( (ac+k-1)*T+1:anni*T , : )];
        MV(:,:,k)=[m( 1:(k-1)*T , : ) ; m( (ac+k-1)*T+1:anni*T , : )];
        SDV(:,:,k)=[sd( 1:(k-1)*T , : ) ; sd( (ac+k-1)*T+1:anni*T , : )];
        
    end
    
end

%% ANN in cross-validazione per trovare i parametri ottimali

%Utilizziamo la crossvalidazione per capire il numero ottimale di neuroni

n_max_neu = 5; %numero massimo di neuroni ###### da cambiare
R2_opt_c = nan*zeros(n_max_neu,K); 
R2_opt_v = nan*zeros(n_max_neu,K);
N_runs = 5; % ####### da cambiare
    
for k=1:K %per ogni dataset della crossvalidazione
    
    dc=DC(:,:,k);
    xc=XC(:,:,k);
    mc=MC(:,:,k);
    sdc=SDC(:,:,k);
    dv=DV(:,:,k);
    xv=XV(:,:,k);
    mv=MV(:,:,k);
    sdv=SDV(:,:,k);

    X = xc;
    Xv = xv;
    R2_ann_cal = nan*ones(n_max_neu,N_runs);
    R2_ann_val = nan*ones(n_max_neu,N_runs);
    %Il generico elemento (i,j) delle tabelle R2_ann_cal e R2_ann_val 
    %conterrà l'R2 corrispondente a una ANN con i neuroni e j iterazioni
    
    for i=1:n_max_neu %i = current number of neurons
        
        for j = 1:N_runs
            
            %calibration
            net = feedforwardnet(i);
            net = train(net,X(:,1:end-1)',X(:,end)');
            Y = [net(X(:,1:end-1)')'];
            d_ann = Y.*sdc(:,end) + mc(:,end);

            R2_ann_cal(i,j) = 1 - sum((dc(:,end) - d_ann).^2) / sum((dc(:,end)-mc(:,end)).^2);

            %validation
            Yv = [net(Xv(:,1:end-1)')'];
            d_ann_v = Yv.*sdv(:,end) + mv(:,end);

            R2_ann_val(i,j) = 1 - sum((dv(:,end) - d_ann_v).^2) / sum((dv(:,end)-mv(:,end)).^2);
        end
        
    end
    
    for i = 1:n_max_neu
        R2_opt_c(i,k) = max( R2_ann_cal(i,:));
        %R2_opt_c conterrà nella riga i, l'R2 massimo (tra tutte le N_runs)
        %in calibrazione con i-neuroni. La colonna k corrisponde al dataset k-esimo.
    end

    for i = 1:n_max_neu
        R2_opt_v(i,k) = max( R2_ann_val(i,:));
        %R2_opt_v conterrà nella riga i, l'R2 massimo (tra tutte le N_runs)
        %in validazione con i-neuroni. La colonna k corrisponde al dataset k-esimo.
    end
    
a_che_punto_siamo = k/K*100

end

R2_ann_c = mean(R2_opt_c,2) 
%Nel vettore R2_ann_c ogni riga i corrisponde all'R2 medio (in calibrazione) in
%crossvalidazione, con i-neuroni
R2_ann_v = mean(R2_opt_v,2)
%Nel vettore R2_ann_v ogni riga i corrisponde all'R2 medio (in validazione) in
%crossvalidazione, con i-neuroni

[ R2_max_ann , opt_neu] = max(R2_ann_v);
% R2_max_ann = R2 max (con numero ottimale di neuroni)
% opt_neu = numero ottimale di neuroni

%Visualizzo gli R2 con un barplot
barplot=bar(R2_ann_v);
title('Mean R2 in cross-validation for increasing no. of neurons', 'fontsize', 15);
xlabel('no. of neurons', 'fontsize', 13);
ylabel('Mean R2 in cross-validation', 'fontsize', 13)
barplot.FaceColor='flat';
colorData = barplot.CData;
barplot.CData(opt_neu,:) = [1 0 0]; 

%% Trovo gli ordini ottimali con il dataset completo diviso 5-2 

k=1;
    dc=DC(:,:,k);
    xc=XC(:,:,k);
    mc=MC(:,:,k);
    sdc=SDC(:,:,k);
    dv=DV(:,:,k);
    xv=XV(:,:,k);
    mv=MV(:,:,k);
    sdv=SDV(:,:,k);
    
X = xc; 
Xv = xv;
N_runs = 5;
max_lag = 7; %##### da cambiare
num_input = size(dc,2)-1; %(escludo l'ultima colonna che è l'output da stimare)
Ann_ord_opt = nan(max_lag,num_input);
R2_ann_cal = nan*ones(max_lag,N_runs);
R2_ann_val = nan*ones(max_lag,N_runs);

 for l = 1:num_input
     
     X=xc(max_lag:end,l); 
     Xv=xv(max_lag:end,l);
     
     for i = 1:max_lag
         
         if i>1 
         X=[X xc(max_lag-i+1:end-i+1,l)];
         Xv=[Xv xv(max_lag-i+1:end-i+1,l)];
         end
         
        for j = 1:N_runs
            
            %calibration
            
            net = feedforwardnet(opt_neu);
            net = train(net,X',xc(max_lag:end,end)');
            Y = [net(X')'];
            d_ann = Y.*sdc(max_lag:end,end) + mc(max_lag:end,end);

            R2_ann_cal(i,j) = 1 - sum((dc(max_lag:end,end) - d_ann).^2) / sum((dc(max_lag:end,end)-mc(max_lag:end,end)).^2);

            %validation
            
            Yv = [net(Xv')'];
            d_ann_v = Yv.*sdv(max_lag:end,end) + mv(max_lag:end,end);

            R2_ann_val(i,j) = 1 - sum((dv(max_lag:end,end) - d_ann_v).^2) / sum((dv(max_lag:end,end)-mv(max_lag:end,end)).^2);  
        end
        
        Ann_ord_opt(i,l) = max( R2_ann_val(i,:));
        
     end
    
     a_che_punto_siamo = l/num_input*100
     
 end
 
 [r2,ordini] = max(Ann_ord_opt);
 
%% ANN finale con il dataset con gli ordini ottimali divisi 5-2

% costruisco il dataset
flag=0;
for i=1:length(r2)
   for j=1:ordini(i)
        if flag==0
             X=xc(max(ordini)-j+1:end-j+1,i);
             Xv=xv(max(ordini)-j+1:end-j+1,i);
             flag=1;
        end
        if flag==1
             X=[X xc(max(ordini)-j+1:end-j+1,i)];
             Xv=[Xv xv(max(ordini)-j+1:end-j+1,i)];
        end
    end
end
 
 % faccio l'ultima rete neurale

n_max_neu = opt_neu;
N_runs = 1;
R2_ann_cal = zeros(1,N_runs);
R2_ann_val = zeros(1,N_runs);

for j = 1:N_runs
    %calibration
    net = feedforwardnet(n_max_neu);
    net = train(net,X(:,1:end-1)',xc(max(ordini):end,end)');
    Y = [net(X(:,1:end-1)')'];
    d_ann = Y.*sdc(max(ordini):end,end) + mc(max(ordini):end,end);

    R2_ann_cal(j) = 1 - sum((dc(max(ordini):end,end) - d_ann).^2) / sum((dc(max(ordini):end,end)-mc(max(ordini):end,end)).^2);

    %validation
    Yv = [net(Xv(:,1:end-1)')'];
    d_ann_v = Yv.*sdv(max(ordini):end,end) + mv(max(ordini):end,end);

    R2_ann_val(j) = 1 - sum((dv(max(ordini):end,end) - d_ann_v).^2) / sum((dv(max(ordini):end,end)-mv(max(ordini):end,end)).^2);
    
    if R2_ann_val(j) >= max(R2_ann_val)
        net_opt = net;
    end
end

R2_opt = max(R2_ann_val)

%% Analisi errore

Yc = (net_opt(X(:,1:end-1)')');
Yv = (net_opt(Xv(:,1:end-1)')');

err_cal = Yc-X(:,end)
err_val = Yv-Xv(:,end)

hist(err_cal)
hist(err_val)

correlogram(err_val,err_val,10)