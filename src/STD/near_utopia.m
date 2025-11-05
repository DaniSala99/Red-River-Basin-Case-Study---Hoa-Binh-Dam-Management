function indice = near_utopia(Pareto , Hydro1, Hydro2, Floods1, Floods2)
n=size(Pareto,1);
dist = nan*zeros(n,1);
min_dist = sqrt(((Pareto(1,end-3) - Hydro1)/(Hydro1-Hydro2))^2 + ((Pareto(1,end-2) - Floods1)/(Floods1-Floods2))^2);%assegno alla minima distanza un valore sicuramente maggiore alla distanza minima)
indice=1;
for i = 2:n 
    dist(i) = sqrt(((Pareto(i,end-3) - Hydro1)/(Hydro1-Hydro2))^2 + ((Pareto(i,end-2) - Floods1)/(Floods1-Floods2))^2);
    
    if dist(i) < min_dist
        min_dist = dist(i);
        indice = i;
    end
end
end