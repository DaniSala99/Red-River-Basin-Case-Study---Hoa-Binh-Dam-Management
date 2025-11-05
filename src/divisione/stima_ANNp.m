function [stima, r2_top] = stima_ANNp(xc,xv,mc,mv,sdc,sdv,qc,qv,n,N_runs,p)

Xc = xc(p:end-1);
Yc = xc(p+1:end);
Xv = xv(p:end-1);
Yv = xv(p+1:end);
val = zeros(N_runs,p);
r2 = nan*zeros(2,1);
for i=1:N_runs
        %calibrazione
        net = feedforwardnet(n);
        net = train(net,Xc',Yc');
        Y = [Yc(1:p); net(Xc')'];
        d_ann = Y.*sdc + mc;
        %validazione
        Yv = [Yc(1:p); net(Xv')'];
        d_ann_v = Yv.*sdv + mv;
        val(i,1) = 1 - sum((qv - d_ann_v).^2) / sum((qv-mv).^2);
        if val(i,1) >= max(val(:,1))
        stima_temp = [ d_ann; d_ann_v ];
        r2(1) = val(i,1);
        end
end
stima = stima_temp;
r2(2) = r2(1);

for i = 2:p
Xc = [xc(p-i+1:end-i) Xc];
Xv = [xv(p-i+1:end-i) Xv];
for j=1:N_runs
        %calibrazione
        net = feedforwardnet(n);
        net = train(net,Xc',Yc');
        Y = [Yc(1:p); net(Xc')'];
        d_ann = Y.*sdc + mc;
        %validazione
        Yv = [Yv(1:p); net(Xv')'];
        d_ann_v = Yv.*sdv + mv;
        val(j,i) = 1 - sum((qv - d_ann_v).^2) / sum((qv-mv).^2);
   if val(j,i) >= max(val(:,i))
       stima_temp = [ d_ann; d_ann_v ];
       r2(1) = val(j,i); 
   end
end
if r2(1) >= r2(2)
    stima=stima_temp;
    r2(2)=r2(1);
end
end
r2_top=r2(2);
end
