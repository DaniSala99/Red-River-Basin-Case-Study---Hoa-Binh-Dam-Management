# Quick Start Guide

Get up and running with the Red River Basin Dam Management project in 10 minutes!

## ⚡ Quick Installation

### 1. Prerequisites Check
```matlab
% Verify MATLAB version (need R2020a+)
version

% Verify required toolboxes
ver
```

Required toolboxes:
- ✅ Optimization Toolbox
- ✅ Deep Learning Toolbox (Neural Network)
- ✅ Statistics and Machine Learning Toolbox

### 2. Clone Repository
```bash
git clone https://github.com/yourusername/red-river-basin-dam-management.git
cd red-river-basin-dam-management
```

### 3. Setup MATLAB Paths
```matlab
% In MATLAB, navigate to project root, then run:
addpath(genpath(pwd))  % Add all subfolders to path
savepath  % Save for future sessions (optional)
```

---

## 🎮 Running Your First Analysis

### Option 1: Forecast Model (Quick - 5 minutes)

```matlab
% Navigate to src/main
cd src/main

% Open and run the ANN forecast script
open FILE_CONCEGNA_FINALE.m

% Run section by section using Ctrl+Enter (Cmd+Enter on Mac)
% Or run all: press F5
```

**Expected output**:
- Bar plot showing R² for different neuron numbers
- Time series plot comparing real vs predicted inflow
- Error distribution histogram
- Correlogram of residuals

**Interpret results**:
- Optimal neuron number is highlighted in red
- R² ~ 0.63 is good performance
- Error should be approximately Gaussian

### Option 2: Policy Optimization (Longer - 30 minutes)

```matlab
% Navigate to src/main
cd src/main

% Open the main analysis script
open NRM_Lab_2022_HoaBinh_Script.m

% Run MANDATORY PART sections only (faster):
% 1. Data import and model setup (Ctrl+Enter)
% 2. Simulation of STD policy (Ctrl+Enter)
% 3. Optimization via NSGAII of STD (Ctrl+Enter) - WAIT ~20 min
% 4. Interesting solutions of STD (Ctrl+Enter)
```

**Expected output**:
- Standard Operating Policy curve
- Time series of reservoir level, release, and Hanoi flooding
- Initial vs final population scatter plot
- Pareto front with interesting solutions marked

**Interpret results**:
- Blue points = initial random policies
- Red points = optimized policies
- Bottom-left is better (less flooding, more hydropower)

---

## 📊 Example Workflow

### Scenario: "I want to compare different dam management strategies"

```matlab
%% 1. Load data and setup
clear; close all; clc;
addpath ./data ./sim ./STD ./RBF ./NSGA2

qq = model_setup([1 11 1994], [31 12 2005]);
h0 = 104;  % Initial water level (m)

%% 2. Test current policy
theta_sop = [89, 107, 1200, 2500, 5000];  % SOP parameters
[h, u, r, ht_HN] = simHB(qq, h0, theta_sop);

% Calculate performance
gt_hyd = g_hydropower(r(sys_param.warmup+1:end), h(sys_param.warmup+1:end));
gt_flo = g_flood(ht_HN(sys_param.warmup:end)*100);
J_hyd = mean(gt_hyd)  % Hydropower objective
J_flo = mean(gt_flo)  % Flooding objective

%% 3. Optimize policy
global opt_inputs;
opt_inputs.qq = qq;
opt_inputs.h0 = h0;
opt_inputs.theta = theta_sop;

[chr0, chrF] = nsga_2(70, 50, 2, 5, ...
                      [84 102 0 0 100], ...
                      [101 117 20000 20000 10000]);

%% 4. Extract best solutions
Pareto_fin = find_pareto_fin(chrF);
[~, idx_hydro] = min(Pareto_fin(:,end-3));  % Best hydropower
[~, idx_flood] = min(Pareto_fin(:,end-2));  % Best flood control

fprintf('Best Hydropower: J_hyd = %.2e, J_flo = %.2f\n', ...
        -Pareto_fin(idx_hydro,end-3), Pareto_fin(idx_hydro,end-2));
fprintf('Best Flood Control: J_hyd = %.2e, J_flo = %.2f\n', ...
        -Pareto_fin(idx_flood,end-3), Pareto_fin(idx_flood,end-2));
```

---

## 🐛 Troubleshooting

### Problem: "Undefined function or variable"

**Likely cause**: Missing subfolders in path

**Solution**:
```matlab
addpath ./data ./sim ./STD ./RBF ./NSGA2 ./divisione
% Or recursively add all:
addpath(genpath(pwd))
```

### Problem: "Out of memory"

**Likely cause**: Large population size or long time series

**Solution**:
```matlab
% Reduce population/generations
pop = 20;  % instead of 70
gen = 20;  % instead of 50

% Clear variables between runs
clear; close all; clc;
```

### Problem: "Optimization takes too long"

**Solutions**:
1. **Reduce problem size**:
   ```matlab
   pop = 30;  % Smaller population
   gen = 20;  % Fewer generations
   ```

2. **Use parallel computing** (if toolbox available):
   ```matlab
   % NSGA-II can be parallelized - check algorithm settings
   ```

3. **Run overnight**: Some analyses are meant to run for hours

### Problem: Neural network won't train

**Solutions**:
1. Check data format:
   ```matlab
   size(X)  % Should be [N_samples, N_features]
   any(isnan(X(:)))  % Should be false
   ```

2. Reduce network complexity:
   ```matlab
   n_neurons = 2;  % Instead of 3 or more
   ```

3. Check training data:
   ```matlab
   plot(X(:,1))  % Visualize input data
   plot(X(:,end))  % Visualize output data
   ```

---

## 📈 Next Steps

After your first successful run:

1. **📖 Read the full README**: Understand theory and methodology
2. **📄 Check the report**: Detailed analysis in `docs/Report_group26.pdf`
3. **🔧 Modify parameters**: Experiment with different settings
4. **📊 Analyze results**: Use MATLAB's plotting tools
5. **🤝 Contribute**: See `CONTRIBUTING.md` for guidelines

---

## 💡 Tips for Best Results

### Performance
- Start with small test runs (fewer generations)
- Use section execution (Ctrl+Enter) to run parts independently
- Save workspace after long runs: `save('results.mat')`

### Analysis
- Always visualize results: `figure; plot(...)`
- Compare multiple solutions on Pareto front
- Check residuals for model validation
- Look for patterns in time series plots

### Reproducibility
- Set random seed for consistency:
  ```matlab
  rng(42);  % Fixed random seed
  ```
- Document your parameter choices
- Save important results: `saveas(gcf, 'figure.png')`

---

## 🎯 Quick Commands Reference

```matlab
% Path management
addpath(genpath(pwd))  % Add all subfolders
savepath               % Save path for future

% Running scripts
run('script_name.m')   % Run entire script
% Ctrl+Enter           % Run current section

% Data management
load('data.mat')       % Load data
save('results.mat')    % Save workspace

% Plotting
figure;                % New figure window
hold on;               % Keep previous plot
legend('label');       % Add legend
saveas(gcf, 'plot.png')% Save figure

% Help
help function_name     % Function documentation
doc function_name      % Detailed documentation
which function_name    % Find function location
```

---

## 📬 Need Help?

- 📝 Open an issue on GitHub
- 📧 Contact maintainers
- 📚 Check documentation in `docs/`
- 💬 Join discussions

---

**You're all set! Happy analyzing! 🚀**
