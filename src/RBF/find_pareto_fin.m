function P = find_pareto_fin(chr)
[m,k]=min(chr(:,end-2))
P = chr(k,:);

for i=1:size(chr,1)
    in=0;
    for j=1:size(chr,1)
    if chr(i,end-3)>chr(j,end-3) && chr(i,end-2)>chr(j,end-2)
        in=1;
    end
    end
    if in==0
        P =[P;chr(i,:)];
    end
end

end