%% MANDATORY PART
% Run once a time each section from 'MANDATORY PART' to 'Interesting solution of STD' 
clear
close all
clc

%% Data import and model setup
addpath ./data
addpath ./sim
addpath ./STD
addpath ./divisione
global sys_param; % global structure containing model parameters

qq = model_setup([1 11 1994], [31 12 2005]); % structure containing trajectories of Hoa Binh inflow and flow of Thao and Lo rivers

%% simulation of STD policy

%SYS_PARAM

%lsv --> level, surface, volume(storage)
%max_rel --> h(level), max release(flow)


% set policy parameters
theta = [89,107,1200,2500,5000]; % h1, h2, m1, m2, w

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

figure;
plot(hh,rr, 'LineWidth', 2);
hold on
pl1=plot(hh, rmin, '--', 'Color', [1 0 0],'LineWidth', 1);
hold on
plot(hh, rmax, '--', 'Color', [0 1 0],'LineWidth', 1);
xlabel('Level [m]','Fontsize', 13);
ylabel('Release [m3/s]', 'Fontsize', 13);
title("Standard Operating Policy", 'Fontsize', 16);
leg=legend('Policy','Min release','Max release', 'Fontsize', 15);

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qq, h0,theta);


% visualize trajectories
date_day = datetime(1995,1,1):datetime(2005,12,31);
figure; 
subplot(221); plot(date_day, qq.q_Da(sys_param.warmup:end),'LineWidth', 1)% inflow
title('Da streamflow', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Streamflow [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
subplot(222);plot(date_day, h(sys_param.warmup+1:end),'LineWidth', 1)
title('Level of reservoir', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
subplot(223); plot(date_day, r(sys_param.warmup+1:end),'LineWidth', 1)
title('Release', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
subplot(224); plot(date_day, ht_HN(sys_param.warmup+1:end),'LineWidth', 1)
yline(9.5,"r-.",'LineWidth', 1);
title('Level in Hanoi', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level[m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 15)

% compute policy perfomance
gt_hyd = g_hydropower( r(sys_param.warmup+1:end), h(sys_param.warmup+1:end) );
gt_flo = g_flood( ht_HN(sys_param.warmup:end)*100 ) ;
Jhyd_std = mean(gt_hyd);
Jflo_std = mean(gt_flo);

%% Optimization via NSGAII of STD

addpath ./NSGA2
global opt_inputs ;
opt_inputs.qq = qq;
opt_inputs.h0 = h0 ;
opt_inputs.theta= theta;
pop = 70; % number of individuals in the population
gen = 50; % number of generations
M = 2; % number of objectives
V = 5; % number of decision variables (policy parameters)
min_range = [84 102 0 0 100]; % minimum value of each parameter
max_range = [101 117 20000 20000 10000 ]; % maximum value of each parameter
[ chr0, chrF ] = nsga_2(pop,gen,M,V,min_range,max_range) ;

%find pareto front
Pareto_in = find_pareto_in(chr0);
Pareto_fin = find_pareto_fin(chrF);

figure
plot( chr0(:,end), chr0(:,end-1), 'o', 'MarkerFaceColor', 'b');
hold on
plot( chrF(:,end-2), chrF(:,end-3), 'ro', 'MarkerFaceColor', 'r');
hold on
plot(Jflo_std,-Jhyd_std,'d','Color', [0 0 0],'LineWidth', 2)
legend('Initial population', 'Final population','SOP', 'Fontsize', 20);
xlim([0 2000]);
tit = title('Initial and final population - nsga2 on Standard Operating Policy', 'Fontsize', 17)
xlab = xlabel('J flooding', 'FontWeight', 'bold', 'Fontsize', 13)
ylabel('J hydropower', 'FontWeight', 'bold', 'Fontsize', 13)

%% Interesting solutions of STD
% select few "interesting" solutions (best hydropower, best flooding, one
% compromise) and visualize the policy along with the simulated
% trajectorie
ob_sol_std= nan*zeros(3,3);

%HYDRO-MAX
[Hydro_max, indice] = min( Pareto_fin(:,end-3));
Hydro_min = max( Pareto_fin(:,end-3));
theta = build_theta(Pareto_fin , indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_hydro = rr;

% set initial condition and run simulation

h0=104; %m
[h_hy, u_hy, r_hy,ht_HN_hy] = simHB(qq, h0,theta);

% visualize trajectories
ob_sol_std(1,1)= "Max Hydro";
ob_sol_std(1,2)= Pareto_fin(indice,end-2);
ob_sol_std(1,3)= Pareto_fin(indice,end-3);

%FLOODS-MIN
[Flood_min, indice] = min( Pareto_fin(:,end-2));
Flood_max = max( Pareto_fin(:,end-2));
theta = build_theta(Pareto_fin , indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_floods = rr;
% set initial condition and run simulation

h0=104; %m
[h_fl, u_fl, r_fl,ht_HN_fl] = simHB(qq, h0,theta);

ob_sol_std(2,1)= "Min Floods";
ob_sol_std(2,2)= Pareto_fin(indice,end-2);
ob_sol_std(2,3)= Pareto_fin(indice,end-3);

% UTOPIA NEAREST

indice = near_utopia(Pareto_fin, Hydro_max, Hydro_min, Flood_min, Flood_max);
theta = build_theta(Pareto_fin, indice);


% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_utopia = rr;

% set initial condition and run simulation

h0=104; %m
[h_ut, u_ut, r_ut,ht_HN_ut] = simHB(qq, h0,theta);

ob_sol_std(3,1)= "Utopia nearest";
ob_sol_std(3,2)= Pareto_fin(indice,end-2);
ob_sol_std(3,3)= Pareto_fin(indice,end-3);

%plot Pareto di confronto
x_fin_std=sortrows(Pareto_fin, size(Pareto_fin,2)-3);
x_in_std=sortrows(Pareto_in , size(Pareto_in,2)-1);
figure
plot(x_fin_std(:,end-2),x_fin_std(:,end-3),'-.+g','LineWidth',2)
hold on
plot(x_in_std(:,end),x_in_std(:,end-1),'-.+b','LineWidth',2)
hold on
plot(ob_sol_std(1,2),ob_sol_std(1,3),'ro','LineWidth',2)
hold oN
plot(ob_sol_std(2,2),ob_sol_std(2,3),'bo','LineWidth',2)
hold on
plot(ob_sol_std(3,2),ob_sol_std(3,3),'go','LineWidth',2)
xlab = xlabel('J flooding', 'FontWeight', 'bold', 'Fontsize', 13)
ylabel('J hydropower', 'FontWeight', 'bold', 'Fontsize', 13)
legend("Final Pareto","Initial Pareto","Best Hydropower","Best Floods","Best Compromise",'Fontsize', 15)
title("Pareto front",'Fontsize', 16)

figure
plot(hh,rr_floods, 'LineWidth',2);
hold on
plot(hh,rr_hydro, 'LineWidth',2);
hold on
plot(hh,rr_utopia, 'LineWidth',2);
hold on
plot(hh, rmin, '--', 'Color', [0 0 0],'LineWidth',1);
hold on
plot(hh, rmax, '--', 'Color', [0 0 0],'LineWidth',1);
legend("min flooding","max hydropower","best compromise",'Fontsize', 15)
title('Standard Operating Policy - 3 optimal solutions', 'Fontsize', 16)
xlabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m3/s]','FontWeight', 'bold', 'Fontsize', 13)

% LEVEL IN HANOI (ht_HN)

date_day = datetime(1995,1,1):datetime(2005,12,31);
figure; 
subplot(221); plot(date_day, ht_HN(sys_param.warmup+1:end),'LineWidth', 1)
yline(9.5,"r-.");
title('Level in Hanoi - SOP', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(222);plot(date_day, ht_HN_fl(sys_param.warmup+1:end),'LineWidth', 1)
yline(9.5,"r-.");
title('Level in Hanoi - Min Flooding', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(223);plot(date_day, ht_HN_hy(sys_param.warmup+1:end),'LineWidth', 1)
yline(9.5,"r-.");
title('Level in Hanoi - Max Hydropower', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(224);plot(date_day, ht_HN_ut(sys_param.warmup+1:end),'LineWidth', 1)
yline(9.5,"r-.");
title('Level in Hanoi - Best compromise', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level[m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])

% H RESERVOIR 
figure; 
subplot(221); plot(date_day, h(sys_param.warmup+1:end),'LineWidth', 1)
title('Level of reservoir - SOP', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(222);plot(date_day, h_fl(sys_param.warmup+1:end),'LineWidth', 1)
title('Level of reservoir - Min Flooding', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(223);plot(date_day, h_hy(sys_param.warmup+1:end),'LineWidth', 1)
title('Level of reservoir - Max Hydropower', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(224);plot(date_day, h_ut(sys_param.warmup+1:end),'LineWidth', 1)
title('Level of reservoir - Best compromise', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])

% RELEASE 
figure; 
subplot(221); plot(date_day, r(sys_param.warmup+1:end),'LineWidth', 1)
title('Release - SOP', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(222);plot(date_day, r_fl(sys_param.warmup+1:end),'LineWidth', 1)
title('Release - Min Flooding', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(223);plot(date_day, r_hy(sys_param.warmup+1:end),'LineWidth', 1)
title('Release - Max Hydropower', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(224);plot(date_day, r_ut(sys_param.warmup+1:end),'LineWidth', 1)
title('Release - Best compromise', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])

%% NON MANDATORY PART
clear
close all
clc

%% Data import and model setup
addpath ./data
addpath ./sim
addpath ./STD
addpath ./divisione
global sys_param; % global structure containing model parameters

qq = model_setup([1 11 1994], [31 12 2005]); % structure containing trajectories of Hoa Binh inflow and flow of Thao and Lo rivers

date_day = datetime(1994,11,1):datetime(2005,12,31);
date_cal = find(date_day == datetime(2001,12,31));

qqc.q_Da = qq.q_Da(1:date_cal);
qqc.q_YB = qq.q_YB(1:date_cal);
qqc.q_VQ = qq.q_VQ(1:date_cal);

qqv.q_Da = qq.q_Da(date_cal+1:end);
qqv.q_YB = qq.q_YB(date_cal+1:end);
qqv.q_VQ = qq.q_VQ(date_cal+1:end);


figure
subplot(221);
xline(datetime(2001,12,31),"r-.",'LineWidth',1);
hold on 
plot(date_day,qq.q_Da,'LineWidth',1)
title('Da steamflow', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Streamflow [m3/s]','FontWeight', 'bold', 'Fontsize', 13)
leg=legend('Calibration limit', 'Fontsize', 12);
ylim([0 18000])
subplot(222);
xline(datetime(2001,12,31),"r-.",'LineWidth',1);
hold on 
plot(date_day,qq.q_YB,'LineWidth',1)
title('Thao steamflow', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Streamflow [m3/s]','FontWeight', 'bold', 'Fontsize', 13)
leg=legend('Calibration limit', 'Fontsize', 12);
ylim([0 18000])
subplot(223);
xline(datetime(2001,12,31),"r-.",'LineWidth',1);
hold on 
plot(date_day,qq.q_VQ,'LineWidth',1)
title('Lo steamflow', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Streamflow [m3/s]','FontWeight', 'bold', 'Fontsize', 13)
leg=legend('Calibration limit', 'Fontsize', 12);
ylim([0 18000])
subplot(224);
xline(datetime(2001,12,31),"r-.",'LineWidth',1);
hold on
plot(date_day,qq.q_VQ+qq.q_YB,'LineWidth',1)
title('Lo + Thao steamflow', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Streamflow [m3/s]','FontWeight', 'bold', 'Fontsize', 13)
leg=legend('Calibration limit', 'Fontsize', 12);
ylim([0 18000])

%% simulation of STD policy

%SYS_PARAM

%lsv --> level, surface, volume(storage)
%max_rel --> h(level), max release(flow)


% set policy parameters
theta = [89,107,1200,2500,5000]; % h1, h2, m1, m2, w

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr_std(i) = min( max(rmin(i), uu(i)), rmax(i) );
end


% plot the BASELINE POLICY
figure;
plot(hh,rr_std, 'LineWidth', 2);
hold on
pl1=plot(hh, rmin, '--', 'Color', [1 0 0],'LineWidth', 1);
hold on
plot(hh, rmax, '--', 'Color', [0 1 0],'LineWidth', 1);
xlabel('Level [m]','Fontsize', 13);
ylabel('Release [m3/s]', 'Fontsize', 13);
title("Standard Operating Policy", 'Fontsize', 16);
leg=legend('Policy','Min release','Max release', 'Fontsize', 10);

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qq, h0,theta);

h_std=h;
u_std=u;
r_std=r;
ht_HN_std=ht_HN;

% visualize trajectories
date_day = datetime(1995,1,1):datetime(2005,12,31);
figure; 
subplot(221); plot(date_day, qq.q_Da(sys_param.warmup:end))% inflow
title("Da stream inflow")
subplot(222); plot(date_day, h(sys_param.warmup+1:end))% HB level
title("Level of reservoir");
subplot(223); plot(date_day, r(sys_param.warmup+1:end))% HB release
title("Release of HoaBin");
subplot(224); plot(date_day, ht_HN(sys_param.warmup+1:end))% level in Hanoi
yline(9.5,"r-.");
title("Level in Hanoi");

% compute policy perfomance
gt_hyd = g_hydropower( r(sys_param.warmup+1:end), h(sys_param.warmup+1:end) );
gt_flo = g_flood( ht_HN(sys_param.warmup:end)*100 ) ;
Jhyd_std = mean(gt_hyd);
Jflo_std = mean(gt_flo);

%% Optimization via NSGAII of STD

addpath ./NSGA2
global opt_inputs ;
opt_inputs.qq = qqc;
opt_inputs.h0 = h0 ;
opt_inputs.theta= theta;
pop = 70; % number of individuals in the population
gen = 50; % number of generations
M = 2; % number of objectives
V = 5; % number of decision variables (policy parameters)
min_range = [84 102 0 0 100]; % minimum value of each parameter
max_range = [101 117 10000 10000 10000 ]; % maximum value of each parameter
[ chr0, chrF ] = nsga_2(pop,gen,M,V,min_range,max_range) ;

% visualize the initial and final population and the performance of the SOP
figure
plot( chr0(:,end), chr0(:,end-1), 'o', 'MarkerFaceColor', 'b');
% col 5 = hydropower, col 6 = flood
hold on
plot( chrF(:,end-2), chrF(:,end-3), 'ro', 'MarkerFaceColor', 'r');
leg=legend('Initial', 'Optimized', 'Fontsize', 20);
xlim([200 1000]);
tit = title('Initial and final population - nsga2 on Standard Operating Policy', 'Fontsize', 17)
xlab = xlabel('J flooding', 'FontWeight', 'bold', 'Fontsize', 13)
ylabel('J hydropower', 'FontWeight', 'bold', 'Fontsize', 13)

popin_Jflo_std=chr0(:,end);
popin_Jhyd_std=chr0(:,end-1);
popfin_Jflo_std=chrF(:,end-2);
popfin_Jhyd_std=chrF(:,end-3);

%find pareto front
Pareto_in = find_pareto_in(chr0);
Pareto_fin = find_pareto_fin(chrF);

Pareto_in_flood_std = Pareto_in(:,end);
Pareto_in_hydro_std = Pareto_in(:,end-1);
Pareto_fin_flood_std = Pareto_fin(:,end-2);
Pareto_fin_hydro_std = Pareto_fin(:,end-3);

%% Simulation in validation of STD

chrV = nan * zeros(size(chrF,1),size(chrF,2));

for i = 1:size(chrF,1)
    
    theta = chrF(i,1:end-4);
    chrV(i,1:end-4) = theta;
    
    h0=104;
    [h, u, r,ht_HN] = simHB(qqv, h0,theta);
    
    gt_hyd = g_hydropower( r(2:end), h(2:end) );
    gt_flo = g_flood( ht_HN(2:end)*100 ) ;
    chrV(i,end-3) = -mean(gt_hyd);
    chrV(i,end-2) = mean(gt_flo);
end
Pareto_val = find_pareto_fin(chrV);

%% Interesting solutions of STD cal
% select few "interesting" solutions (best hydropower, best flooding, one
% compromise) and visualize the policy along with the simulated
% trajectorie
ob_sol_std_cal= nan*zeros(3,3);

%HYDRO-MAX
[Hydro_max, indice] = min( Pareto_fin(:,end-3));
Hydro_min = max( Pareto_fin(:,end-3));
theta = build_theta(Pareto_fin , indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_hydro = rr;
rr_std_hydro = rr;

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqc, h0,theta);

h_std_hydro_cal=h;
u_std_hydro_cal=u;
r_std_hydro_cal=r;
ht_HN_std_hydro_cal=ht_HN;

ob_sol_std_cal(1,1)= "Max Hydro";
ob_sol_std_cal(1,2)= Pareto_fin(indice,end-2);
ob_sol_std_cal(1,3)= Pareto_fin(indice,end-3);

%FLOODS-MIN

[Flood_min, indice] = min( Pareto_fin(:,end-2));
Flood_max = max( Pareto_fin(:,end-2));
theta = build_theta(Pareto_fin , indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_floods = rr;
rr_std_floods = rr;

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqc, h0,theta);

h_std_flood_cal=h;
u_std_flood_cal=u;
r_std_flood_cal=r;
ht_HN_std_flood_cal=ht_HN;

ob_sol_std_cal(2,1)= "Min Floods";
ob_sol_std_cal(2,2)= Pareto_fin(indice,end-2);
ob_sol_std_cal(2,3)= Pareto_fin(indice,end-3);

% UTOPIA NEAREST

indice = near_utopia(Pareto_fin, Hydro_max, Hydro_min, Flood_min, Flood_max);
theta = build_theta(Pareto_fin, indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_utopia = rr;
rr_std_utopia = rr;

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqc, h0,theta);

h_std_utopia_cal=h;
u_std_utopia_cal=u;
r_std_utopia_cal=r;
ht_HN_std_utopia_cal=ht_HN;

ob_sol_std_cal(3,1)= "Utopia nearest";
ob_sol_std_cal(3,2)= Pareto_fin(indice,end-2);
ob_sol_std_cal(3,3)= Pareto_fin(indice,end-3);

% LEVEL IN HANOI (ht_HN)

date_day = datetime(1995,1,1):datetime(2001,12,31);
figure; 
subplot(221); plot(date_day, ht_HN_std(sys_param.warmup+1:date_cal+1))
yline(9.5,"r-.");
title('Level in Hanoi - SOP Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(222);plot(date_day, ht_HN_std_flood_cal(sys_param.warmup+1:end))
yline(9.5,"r-.");
title('Level in Hanoi - Min Flooding Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(223);plot(date_day, ht_HN_std_hydro_cal(sys_param.warmup+1:end))
yline(9.5,"r-.");
title('Level in Hanoi - Max Hydropower Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(224);plot(date_day, ht_HN_std_utopia_cal(sys_param.warmup+1:end))
yline(9.5,"r-.");
title('Level in Hanoi - Best compromise Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level[m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])

% H RESERVOIR 
figure; 
subplot(221); plot(date_day, h_std(sys_param.warmup+1:date_cal+1))
title('Level of reservoir - SOP Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(222);plot(date_day, h_std_flood_cal(sys_param.warmup+1:end))
title('Level of reservoir - Min Flooding Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(223);plot(date_day, h_std_hydro_cal(sys_param.warmup+1:end))
title('Level of reservoir - Max Hydropower Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(224);plot(date_day, h_std_utopia_cal(sys_param.warmup+1:end))
title('Level of reservoir - Best compromise Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])

% RELEASE 
figure; 
subplot(221); plot(date_day, r_std(sys_param.warmup+1:date_cal+1))
title('Release - SOP Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(222);plot(date_day, r_std_flood_cal(sys_param.warmup+1:end))
title('Release - Min Flooding Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(223);plot(date_day, r_std_hydro_cal(sys_param.warmup+1:end))
title('Release - Max Hydropower Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(224);plot(date_day, r_std_utopia_cal(sys_param.warmup+1:end))
title('Release - Best compromise Calibration', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])

figure
plot(hh,rr_std_floods, 'LineWidth',2);
hold on
plot(hh,rr_std_hydro, 'LineWidth',2);
hold on
plot(hh,rr_std_utopia, 'LineWidth',2);
hold on
pl1=plot(hh, rmin, '--', 'Color', [0 0 0]);
hold on
plot(hh, rmax, '-.', 'Color', [0 0 0]);
legend("min flooding","max hydropower","nearest utopia")
title('Standard Operating Policy - 3 optimal solutions', 'Fontsize', 16)
xlabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m3/s]','FontWeight', 'bold', 'Fontsize', 13)

%% Interesting solutions of STD val
% select few "interesting" solutions (best hydropower, best flooding, one
% compromise) and visualize the policy along with the simulated
% trajectorie
ob_sol_std_val= nan*zeros(3,3);

%HYDRO-MAX
[Hydro_max, indice] = min( Pareto_val(:,end-3));
Hydro_min = max( Pareto_val(:,end-3));
theta = build_theta(Pareto_val , indice);

% set initial condition and run simulation

h0=104;
[h, u, r,ht_HN] = simHB(qqv, h0,theta);

h_std_hydro_val=h;
u_std_hydro_val=u;
r_std_hydro_val=r;
ht_HN_std_hydro_val=ht_HN;

ob_sol_std_val(1,1)= "Max Hydro";
ob_sol_std_val(1,2)= Pareto_val(indice,end-2);
ob_sol_std_val(1,3)= Pareto_val(indice,end-3);

%FLOODS-MIN

[Flood_min, indice] = min( Pareto_val(:,end-2));
Flood_max = max( Pareto_val(:,end-2));
theta = build_theta(Pareto_fin , indice);

h0=104; %m
[h, u, r,ht_HN] = simHB(qqv, h0,theta);

h_std_flood_val=h;
u_std_flood_val=u;
r_std_flood_val=r;
ht_HN_std_flood_val=ht_HN;

ob_sol_std_val(2,1)= "Min Floods";
ob_sol_std_val(2,2)= Pareto_val(indice,end-2);
ob_sol_std_val(2,3)= Pareto_val(indice,end-3);

% UTOPIA NEAREST

indice = near_utopia(Pareto_val, Hydro_max, Hydro_min, Flood_min, Flood_max);
theta = build_theta(Pareto_val, indice);

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqv, h0,theta);

h_std_utopia_val=h;
u_std_utopia_val=u;
r_std_utopia_val=r;
ht_HN_std_utopia_val=ht_HN;

ob_sol_std_val(3,1)= "Utopia nearest";
ob_sol_std_val(3,2)= Pareto_val(indice,end-2);
ob_sol_std_val(3,3)= Pareto_val(indice,end-3);

%plot Pareto di confronto
x_fin_std=sortrows(Pareto_fin, size(Pareto_fin,2)-3);
x_in_std=sortrows(Pareto_in , size(Pareto_in,2)-1);
x_fin_std_val=sortrows(Pareto_val, size(Pareto_val,2)-3);
figure
plot(x_fin_std(:,end-2),x_fin_std(:,end-3),'g','LineWidth', 2)
hold on
plot(x_in_std(:,end),x_in_std(:,end-1),'b','LineWidth', 2)
hold on
plot(x_fin_std_val(:,end-2),x_fin_std_val(:,end-3),'r','LineWidth', 2)
hold on
plot(ob_sol_std_cal(1,2),ob_sol_std_cal(1,3),'bo','MarkerFaceColor', 'b')
hold on
plot(ob_sol_std_cal(2,2),ob_sol_std_cal(2,3),'ro', 'MarkerFaceColor', 'r')
hold on
plot(ob_sol_std_cal(3,2),ob_sol_std_cal(3,3),'yo', 'MarkerFaceColor', 'y')
hold on
plot(Jflo_std,-Jhyd_std,'o')
hold on
plot(ob_sol_std_val(1,2),ob_sol_std_val(1,3),'bo', 'MarkerFaceColor', 'b')
hold on
plot(ob_sol_std_val(2,2),ob_sol_std_val(2,3),'ro', 'MarkerFaceColor', 'r')
hold on
plot(ob_sol_std_val(3,2),ob_sol_std_val(3,3),'yo', 'MarkerFaceColor', 'y')
legend("Final Calibration Pareto","Initial Calibration Pareto","Validation Pareto","Best Hydro","Best Floods","Nearest Utopia","SOP")
title("Pareto front", 'Fontsize', 16);
xlabel('J flooding', 'FontWeight', 'bold', 'Fontsize', 13)
ylabel('J hydropower', 'FontWeight', 'bold', 'Fontsize', 13)

% LEVEL IN HANOI (ht_HN)

date_day = datetime(2002,1,1):datetime(2005,12,31);
figure; 
subplot(221); plot(date_day, ht_HN_std(date_cal+2:end))
yline(9.5,"r-.");
title('Level in Hanoi - SOP Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(222);plot(date_day, ht_HN_std_flood_val(2:end))
yline(9.5,"r-.");
title('Level in Hanoi - Min Flooding Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(223);plot(date_day, ht_HN_std_hydro_val(2:end))
yline(9.5,"r-.");
title('Level in Hanoi - Max Hydropower Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])
subplot(224);plot(date_day, ht_HN_std_utopia_val(2:end))
yline(9.5,"r-.");
title('Level in Hanoi - Best compromise Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level[m]','FontWeight', 'bold', 'Fontsize', 13)
legend('Level','Threshold', 'Optimized', 'Fontsize', 10)
ylim([1 17])

% H RESERVOIR 
figure; 
subplot(221); plot(date_day, h_std(date_cal+2:end))
title('Level of reservoir - SOP Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(222);plot(date_day, h_std_flood_val(2:end))
title('Level of reservoir - Min Flooding Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(223);plot(date_day, h_std_hydro_val(2:end))
title('Level of reservoir - Max Hydropower Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])
subplot(224);plot(date_day, h_std_utopia_val(2:end))
title('Level of reservoir - Best compromise Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylim([75 120])

% RELEASE 
figure; 
subplot(221); plot(date_day, r_std(date_cal+2:end))
title('Release - SOP Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(222);plot(date_day, r_std_flood_val(2:end))
title('Release - Min Flooding Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(223);plot(date_day, r_std_hydro_val(2:end))
title('Release - Max Hydropower Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])
subplot(224);plot(date_day, r_std_utopia_val(2:end))
title('Release - Best compromise Validation', 'Fontsize', 16)
xlabel('Time','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m^3/s]','FontWeight', 'bold', 'Fontsize', 13)
ylim([0 20000])

%% Simulation of RBF policy


addpath ./STD -FROZEN
addpath ./RBF

%SYS_PARAM

%lsv --> level, urface, volume(storage)
%max_rel --> h(level), max release(flow)


% set policy parameters
theta = [120,20,117,10,100,5,85,5,15,10,5,2,30000]; % h1, h2, m1, m2, w

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

%plot baseline policy
figure;
plot(hh,rr, 'LineWidth', 2);
hold on
plot(hh,rr_std, 'LineWidth', 2);
hold on
pl1=plot(hh, rmin, '--', 'Color', [0 1 1],'LineWidth', 1);
hold on
plot(hh, rmax, '--', 'Color', [0 1 0],'LineWidth', 1);
xlabel('Level [m]','Fontsize', 13);
ylabel('Release [m3/s]', 'Fontsize', 13);
title("Comparison between SOP and standard RBF", 'Fontsize', 16);
leg=legend('RBF','SOP', 'Min release','Max release', 'Fontsize', 10);


% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qq, h0,theta);

h_rbf=h;
u_rbf=u;
r_rbf=r;
ht_HN_rbf=ht_HN;

% compute policy perfomance
gt_hyd = g_hydropower( r(sys_param.warmup+1:end), h(sys_param.warmup+1:end) );
gt_flo = g_flood( ht_HN(sys_param.warmup:end)*100 ) ;
Jhyd_rbf = mean(gt_hyd);
Jflo_rbf = mean(gt_flo);

%% Optimization via NSGAII of RBF

addpath ./NSGA2
global opt_inputs ;
opt_inputs.qq = qqc;
opt_inputs.h0 = h0 ;
opt_inputs.theta = theta;
pop = 70; % number of individuals in the population
gen = 50; % number of generations
M = 2; % number of objectives
V = 13; % number of decision variables (policy parameters)
min_range = [110 15 90 5 80 2 70 0.5 20 5 1 0.1 10000]; % minimum value of each parameter
max_range = [150 30 120 10 110 5  95 2 30 10 4 2 60000]; % maximum value of each parameter
[ chr0, chrF ] = nsga_2(pop,gen,M,V,min_range,max_range) ;

% visualize the initial and final population and the performance of the RBF
figure
plot( chr0(:,end), chr0(:,end-1), 'o', 'MarkerFaceColor', 'b');
% col 5 = hydropower, col 6 = flood
hold on
plot( chrF(:,end-2), chrF(:,end-3), 'ro', 'MarkerFaceColor', 'r');
leg=legend('Initial', 'Optimized', 'Fontsize', 20);
xlim([0 3000]);
tit=title('Initial and final population - nsga2 on Standard Operating Policy', 'Fontsize', 17)
xlab=xlabel('J flooding', 'FontWeight', 'bold', 'Fontsize', 13)
ylabel('J hydropower', 'FontWeight', 'bold', 'Fontsize', 13)

popin_Jflo_rbf=chr0(:,end);
popin_Jhyd_rbf=chr0(:,end-1);
popfin_Jflo_rbf=chrF(:,end-2);
popfin_Jhyd_rbf=chrF(:,end-3);


%find pareto front
%find pareto front
Pareto_in = find_pareto_in(chr0);
Pareto_fin = find_pareto_fin(chrF);

Pareto_in_flood_rbf = Pareto_in(:,end);
Pareto_in_hydro_rbf = Pareto_in(:,end-1);
Pareto_fin_flood_rbf = Pareto_fin(:,end-2);
Pareto_fin_hydro_rbf = Pareto_fin(:,end-3);

%% Simulation in validation RBF

chrV = nan * zeros(size(chrF,1),size(chrF,2)-2);

for i = 1:size(chrF,1)
    
    theta = chrF(i,1:end-4);
    chrV(i,1:end-2) = theta;
    
    h0=104;
    [h, u, r,ht_HN] = simHB(qqv, h0,theta);
    
    gt_hyd = g_hydropower( r(2:end), h(2:end) );
    gt_flo = g_flood( ht_HN(2:end)*100 ) ;
    chrV(i,end-3) = -mean(gt_hyd);
    chrV(i,end-2) = mean(gt_flo);
end
Pareto_val = find_pareto_fin(chrV);

%% Interesting solutions of RBF cal
% select few "interesting" solutions (best hydropower, best flooding, one
% compromise) and visualize the policy along with the simulated
% trajectorie

ob_sol_rbf_cal= nan*zeros(3,3);

%HYDRO-MAX

[Hydro_max, indice] = min( Pareto_fin(:,end-3));
Hydro_min = max( Pareto_fin(:,end-3));
theta = build_theta(Pareto_fin , indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

%figure;
rr_hydro = rr;
rr_rbf_hydro = rr;

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqc, h0,theta);

h_rbf_hydro_cal=h;
u_rbf_hydro_cal=u;
r_rbf_hydro_cal=r;
ht_HN_rbf_hydro_cal=ht_HN;

ob_sol_rbf_cal(1,1)= "Max Hydro";
ob_sol_rbf_cal(1,2)= Pareto_fin(indice,end-2);
ob_sol_rbf_cal(1,3)= Pareto_fin(indice,end-3);


% FLOODS-MIN

[Flood_min, indice] = min( Pareto_fin(:,end-2));
Flood_max = max( Pareto_fin(:,end-2));
theta = build_theta(Pareto_fin , indice);

% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_floods = rr;
rr_rbf_floods = rr;

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqc, h0,theta);

h_rbf_flood_cal=h;
u_rbf_flood_cal=u;
r_rbf_flood_cal=r;
ht_HN_rbf_flood_cal=ht_HN;

ob_sol_rbf_cal(2,1)= "Min Floods";
ob_sol_rbf_cal(2,2)= Pareto_fin(indice,end-2);
ob_sol_rbf_cal(2,3)= Pareto_fin(indice,end-3);


% UTOPIA NEAREST

indice = near_utopia(Pareto_fin, Hydro_max, Hydro_min, Flood_min, Flood_max);
theta = build_theta(Pareto_fin, indice);


% visualize policy
hh = [75:.5:124] ;
for i=1:length(hh)
    ss(i) = levelToStorage(hh(i));            %storage
    rmin(i) = min_release(ss(i)); %rilascio min in m3/s
    rmax(i) = max_release(ss(i)); %rilascio max in m3/s
    uu(i) = std_operating_policy(hh(i), theta);
    rr(i) = min( max(rmin(i), uu(i)), rmax(i) );
end

rr_utopia = rr;
rr_rbf_utopia = rr;

% set initial condition and run simulation

h0=104; %m
[h, u, r,ht_HN] = simHB(qqc, h0,theta);

h_rbf_utopia_cal=h;
u_rbf_utopia_cal=u;
r_rbf_utopia_cal=r;
ht_HN_rbf_utopia_cal=ht_HN;

ob_sol_rbf_cal(3,1)= "Utopia nearest";
ob_sol_rbf_cal(3,2)= Pareto_fin(indice,end-2);
ob_sol_rbf_cal(3,3)= Pareto_fin(indice,end-3);

%plot Pareto di confronto
x_fin_rbf = sortrows(Pareto_fin, size(Pareto_fin,2)-2);
x_in_rbf = sortrows(Pareto_in , size(Pareto_in,2)-1);
figure
plot(x_fin_rbf(:,end-2),x_fin_rbf(:,end-3),'g')
hold on
plot(x_in_rbf(:,end),x_in_rbf(:,end-1),'b')
hold on
plot(ob_sol_rbf_cal(1,2),ob_sol_rbf_cal(1,3),'ro')
hold on
plot(ob_sol_rbf_cal(2,2),ob_sol_rbf_cal(2,3),'ro')
hold on
plot(ob_sol_rbf_cal(3,2),ob_sol_rbf_cal(3,3),'ro')
legend("Pareto finale","Pareto iniziale","Best Hydro","Best Floods","Best Utopia")
title("Pareto front");

% policy a confronto
fig=figure
pl=plot(hh,rr_floods, 'LineWidth',2);
hold on
plot(hh,rr_hydro, 'LineWidth',2);
hold on
plot(hh,rr_utopia, 'LineWidth',2);
hold on
pl1=plot(hh, rmin, '--', 'Color', [0 0 0]);
hold on
plot(hh, rmax, '-.', 'Color', [0 0 0]);
legend("min flooding","max hydropower","nearest utopia")
title('Radial Basis Function - 3 optimal solutions', 'Fontsize', 16)
xlabel('Level [m]','FontWeight', 'bold', 'Fontsize', 13)
ylabel('Release [m3/s]','FontWeight', 'bold', 'Fontsize', 13)

%% Interesting solutions of RBF val
% select few "interesting" solutions (best hydropower, best flooding, one
% compromise) and visualize the policy along with the simulated
% trajectorie
ob_sol_rbf_val= nan*zeros(3,3);

%HYDRO-MAX
[Hydro_max, indice] = min( Pareto_val(:,end-3));
Hydro_min = max( Pareto_val(:,end-3));
theta = build_theta(Pareto_val , indice);

% set initial condition and run simulation

h0=ht_HN_rbf_hydro_cal(end);
[h, u, r,ht_HN] = simHB(qqv, h0,theta);

h_rbf_hydro_val=h;
u_rbf_hydro_val=u;
r_rbf_hydro_val=r;
ht_HN_rbf_hydro_val=ht_HN;

ob_sol_rbf_val(1,1)= "Max Hydro";
ob_sol_rbf_val(1,2)= Pareto_val(indice,end-2);
ob_sol_rbf_val(1,3)= Pareto_val(indice,end-3);

%FLOODS-MIN

[Flood_min, indice] = min( Pareto_val(:,end-2));
Flood_max = max( Pareto_val(:,end-2));
theta = build_theta(Pareto_fin , indice);

h0=104; %m
[h, u, r,ht_HN] = simHB(qqv, h0,theta);

h_rbf_flood_val=h;
u_rbf_flood_val=u;
r_rbf_flood_val=r;
ht_HN_rbf_flood_val=ht_HN;

ob_sol_rbf_val(2,1)= "Min Floods";
ob_sol_rbf_val(2,2)= Pareto_val(indice,end-2);
ob_sol_rbf_val(2,3)= Pareto_val(indice,end-3);

% UTOPIA NEAREST

indice = near_utopia(Pareto_val, Hydro_max, Hydro_min, Flood_min, Flood_max);
theta = build_theta(Pareto_val, indice);

% set initial condition and run simulation

h0=h_rbf_utopia_cal(end); %m
[h, u, r,ht_HN] = simHB(qqv, h0,theta);

h_rbf_utopia_val=h;
u_rbf_utopia_val=u;
r_rbf_utopia_val=r;
ht_HN_rbf_utopia_val=ht_HN;

ob_sol_rbf_val(3,1)= "Utopia nearest";
ob_sol_rbf_val(3,2)= Pareto_val(indice,end-2);
ob_sol_rbf_val(3,3)= Pareto_val(indice,end-3);

%plot Pareto di confronto
x_fin_rbf=sortrows(Pareto_fin, size(Pareto_fin,2)-3);
x_in_rbf=sortrows(Pareto_in , size(Pareto_in,2)-1);
x_fin_rbf_val=sortrows(Pareto_val, size(Pareto_val,2)-3);
figure
plot(x_fin_rbf(:,end-2),x_fin_rbf(:,end-3),'g','LineWidth', 2)
hold on
plot(x_in_rbf(:,end),x_in_rbf(:,end-1),'b','LineWidth', 2)
hold on
plot(x_fin_rbf_val(:,end-2),x_fin_rbf_val(:,end-3),'r','LineWidth', 2)
hold on
plot(ob_sol_rbf_cal(1,2),ob_sol_rbf_cal(1,3),'bo', 'MarkerFaceColor', 'b')
hold on
plot(ob_sol_rbf_cal(2,2),ob_sol_rbf_cal(2,3),'ro', 'MarkerFaceColor', 'r')
hold on
plot(ob_sol_rbf_cal(3,2),ob_sol_rbf_cal(3,3),'yo', 'MarkerFaceColor', 'y')
hold on
plot(Jflo_rbf,-Jhyd_rbf,'o')
hold on
plot(ob_sol_rbf_val(1,2),ob_sol_rbf_val(1,3),'bo','MarkerFaceColor', 'b')
hold on
plot(ob_sol_rbf_val(2,2),ob_sol_rbf_val(2,3),'ro','MarkerFaceColor', 'r')
hold on
plot(ob_sol_rbf_val(3,2),ob_sol_rbf_val(3,3),'yo','MarkerFaceColor', 'y')
legend("Final Pareto in calibration","Initial Pareto in calibration","Validation Pareto","Best Hydro","Best Floods","Nearest Utopia","SOP")
title("Pareto front")

%% RBF VS STD

figure
plot(x_fin_std_val(:,end-2),x_fin_std_val(:,end-3),'-.og','LineWidth', 2)
hold on
plot(x_fin_rbf_val(:,end-2),x_fin_rbf_val(:,end-3),'-.ob','LineWidth', 2)
%hold on
%plot(x_fin_std(:,end-2),x_fin_std(:,end-3),'-+', 'Color', [0.2 0.6 0.2],'LineWidth', 2)
%hold on
%plot(x_fin_rbf(:,end-2),x_fin_rbf(:,end-3),'-+','Color', [0.3 0.3 0.4],'LineWidth', 2)
hold on
plot(Jflo_rbf,-Jhyd_rbf,'d','Color', [0 0 0],'LineWidth', 2)
hold on
plot(Jflo_std,-Jhyd_std,'p','Color', [0 0 0],'LineWidth', 2)
%hold on
%plot(ob_sol_rbf_val(1,2),ob_sol_rbf_val(1,3),'bo', 'MarkerFaceColor', 'b')
%hold on
%plot(ob_sol_rbf_val(2,2),ob_sol_rbf_val(2,3),'bo', 'MarkerFaceColor', 'b')
%hold on
%plot(ob_sol_rbf_val(3,2),ob_sol_rbf_val(3,3),'bo', 'MarkerFaceColor', 'b')
%hold on
%plot(ob_sol_std_val(1,2),ob_sol_std_val(1,3),'go',  'MarkerFaceColor', 'g')
%hold on
%plot(ob_sol_std_val(2,2),ob_sol_std_val(2,3),'go', 'MarkerFaceColor', 'g')
%hold on
%plot(ob_sol_std_val(3,2),ob_sol_std_val(3,3),'go', 'MarkerFaceColor', 'g')
leg=legend('Piecewise validation', 'RBF validation','Standard radial basis function','SOP', 'Fontsize', 15);
tit = title('Graphical comparison of Pareto fronts', 'Fontsize', 17)
xlab = xlabel('J flooding', 'FontWeight', 'bold', 'Fontsize', 13)
ylabel('J hydropower', 'FontWeight', 'bold', 'Fontsize', 13)

