function [M,y_] = validaz_ARp(u,y,teta)
%dati gli afflussi e i parametri voglio calcolare l'output in validazione
N=length(u);
p=length(teta);
M=zeros(N-p,p);
for i=1:p
  M(:,i) = u(p-i+1:end-i);
end
y_=M*teta;
y_=[y(1:p);y_];
end