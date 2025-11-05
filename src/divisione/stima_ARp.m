function [teta,M,y_] = stima_ARp(x,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
N=length(x);
M=zeros(N-p,p);%variabile da 1 a p giorni prima
y=x(p+1:end);%misure da stimare
for i=1:p
    M(:,i)=x((p-i+1):(end-i));
end
teta=M\y;
y_=M*teta;
y_=[x(1:p);y_];
end
