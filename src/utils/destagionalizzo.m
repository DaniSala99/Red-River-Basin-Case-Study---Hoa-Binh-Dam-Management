function [x,ym,sd] = destagionalizzo(u,T,f)

% Questa funzione destagionalizza una serie temporale

% Input:
% - u = serie temporale
% - T = periodo
% - f = semiampiezza della finestra

% Output:
% - x = serie temporale destagionalizzata
% - ym = media mobile ripetuta per tutta la serie temporale
% - sd = deviazione standard ripetuta per tutta la serie temporale

 [ m , ym ] = moving_average( u , T , f );
 [ s2 , s2_rep ] = moving_average( (u-ym).^2 , T , f );
 sd = sqrt(s2_rep); %deviazione standard ripetuta per tutta la serie temporale
 s = sqrt(s2); %deviazione standard in un solo anno
 
 x = (u-ym)./sd; %destagionalizzo
 
end
