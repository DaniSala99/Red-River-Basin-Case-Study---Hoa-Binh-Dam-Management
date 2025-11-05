function theta = build_theta(Pareto , indice)

n=size(Pareto,2)-2;
theta = [Pareto(indice,1)];

for i=2:n
theta = [theta Pareto(indice,i)]; % h1, h2, m1, m2, w
end
end