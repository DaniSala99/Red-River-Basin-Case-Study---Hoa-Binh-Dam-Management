function P = find_pareto_in(chr)
[m,k]=min(chr(:,end))
P = chr(k,:);

for i=1:size(chr,1)
    in=0;
    for j=1:size(chr,1)
    if chr(i,end-1)>chr(j,end-1) && chr(i,end)>chr(j,end)
        in=1;
    end
    end
    if in==0
        P =[P;chr(i,:)];
    end
end

end