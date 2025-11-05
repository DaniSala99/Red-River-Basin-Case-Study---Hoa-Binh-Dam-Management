# Red River Basin Case Study - Hoa Binh Dam Management

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a+-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

## 📋 Project Overview

This project, developed for the **Natural Resources Management** course (A.Y. 2021/2022) at Politecnico di Milano, addresses the complex multi-objective problem of managing the Hoa Binh dam in Vietnam's Red River Basin. The project aims to optimize dam operations to balance two critical objectives:

- **Hydropower production maximization**
- **Flood control in Hanoi** (capital city with 16 million inhabitants)

### Key Context
- The Red River Basin extends across China, Vietnam, and Laos
- Upstream data from China is unavailable due to strategic reasons
- Hydropower covers ~46% of Vietnam's electrical demand
- Monsoon seasons cause significant flooding risks
- Climate change increases extreme event frequency

---

## 🎯 Project Components

### Part 1: Data-Driven Forecast Model
**Objective**: Forecast the 5-day cumulative inflow to Hoa Binh reservoir

**Key Features**:
- Multiple model structures tested: Linear ARX, ANN, CART, Random Forest
- K-Fold cross-validation (5 years calibration, 2 years validation)
- De-seasonalization based on cyclo-stationary mean/std deviation
- **Best Model**: Artificial Neural Network (3 neurons)
  - R² in calibration: ~0.67
  - R² in validation: ~0.63

**Input Variables**:
- Streamflow measurements: MuongTe, TamDuong, NamGiang, LaiChau, Tabu, HoaBinh
- Precipitation data: MuongTe, TamDuong, Da, BaoLac, BacMe, HaGiang

### Part 2: Policy Optimization via EMODPS
**Objective**: Optimize Hoa Binh dam control policy

**Approach**:
- Evolutionary Multi-Objective Direct Policy Search (EMODPS)
- NSGA-II algorithm (70 individuals, 50 generations)
- Two policy structures compared:
  - Standard Operating Policy (SOP) - 5 parameters
  - Radial Basis Functions (RBF) - 13 parameters

**Objectives**:
- Maximize mean daily hydropower production (kWh/d)
- Minimize mean daily squared water level excess in Hanoi (cm²)

**Results**: RBF policy significantly outperforms SOP, providing better trade-offs between hydropower and flood control.

---

## 📁 Repository Structure

```
red-river-basin-dam-management/
├── README.md                          # This file
├── LICENSE                            # License information
├── docs/
│   └── Report_group26.pdf            # Full technical report
├── src/
│   ├── main/
│   │   ├── NRM_Lab_2022_HoaBinh_Script.m    # Main analysis script
│   │   └── FILE_CONCEGNA_FINALE.m           # ANN forecast model script
│   ├── utils/
│   │   ├── correlogram.m                     # Correlation analysis
│   │   └── destagionalizzo.m                # Deseasonalization function
│   ├── data/                                 # Data files (not included - add your own)
│   ├── sim/                                  # Simulation functions
│   ├── STD/                                  # Standard Operating Policy functions
│   ├── RBF/                                  # Radial Basis Functions
│   ├── NSGA2/                                # NSGA-II optimization algorithm
│   └── divisione/                            # Additional utilities
├── results/                                  # Output figures and results
└── requirements.txt                          # Required MATLAB toolboxes

```

---

## 🚀 Getting Started

### Prerequisites

- MATLAB R2020a or later
- Required MATLAB Toolboxes:
  - Optimization Toolbox
  - Neural Network Toolbox (Deep Learning Toolbox)
  - Statistics and Machine Learning Toolbox

### Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/red-river-basin-dam-management.git
cd red-river-basin-dam-management
```

2. Add necessary MATLAB paths in your script or startup file:
```matlab
addpath('./src/main')
addpath('./src/utils')
addpath('./src/data')
addpath('./src/sim')
addpath('./src/STD')
addpath('./src/RBF')
addpath('./src/NSGA2')
```

### Usage

#### Part 1: Forecast Model (ANN)

Run the ANN forecast model:
```matlab
% Navigate to src/main
cd src/main

% Run the forecast model script
FILE_CONCEGNA_FINALE
```

**Key Steps**:
1. Data loading and deseasonalization
2. K-fold cross-validation partitioning
3. Optimal neuron number identification
4. Optimal lag selection for each input
5. Final model training and validation
6. Error analysis (distribution and autocorrelation)

#### Part 2: Policy Optimization

Run the complete analysis (mandatory + optional parts):
```matlab
% Navigate to src/main
cd src/main

% Run section by section (recommended)
% Use MATLAB's section execution (Ctrl+Enter)
NRM_Lab_2022_HoaBinh_Script

% Or run specific parts:
% - MANDATORY PART: SOP optimization on full dataset (1995-2005)
% - NON MANDATORY PART: Cal/Val split + RBF comparison
```

**Workflow**:
1. **Data Import**: Load streamflow data from Da, Thao, and Lo rivers
2. **SOP Baseline**: Evaluate current Standard Operating Policy
3. **Optimization**: Run NSGA-II on SOP parameters
4. **Pareto Analysis**: Extract and analyze interesting solutions:
   - Best Hydropower
   - Best Flood Control
   - Best Compromise (nearest to utopia)
5. **RBF Comparison** (optional): Compare RBF vs SOP performance

---

## 📊 Key Results

### Forecast Model Performance

| Model | Structure | R² (Cal) | R² (Val) |
|-------|-----------|----------|----------|
| Linear ARX | 8 variables, 24 parameters | ~0.59 | ~0.59 |
| **ANN (Best)** | **3 neurons, 45 inputs** | **~0.67** | **~0.63** |
| CART | Tree-based | ~0.72 | ~0.48 (overfitting) |
| Random Forest | Ensemble | - | ~0.57 |

**Error Characteristics**:
- Mean: ≈0 (unbiased)
- Std: 0.61 (deseasonalized)
- Distribution: Near-Gaussian with slight left skewness (tends to overestimate)
- Autocorrelation: Present up to lag-3 (room for improvement)

### Policy Optimization Results

**Current SOP Performance** (1995-2005):
- Hydropower: 1.69 × 10⁷ kWh/d
- Flooding: 569.6 cm²

**Optimized Solutions**:
- Final Pareto front shows significant improvement over SOP
- RBF policies dominate SOP policies in validation
- Trade-off curve allows decision-makers to choose preferred balance

**Key Insights**:
- RBF's smoothness allows better adaptation to system dynamics
- Cal/Val split essential to avoid overfitting given high inter-annual variability
- Longer time series would improve solution robustness

---

## 🔧 Functions Description

### Core Analysis Scripts

#### `NRM_Lab_2022_HoaBinh_Script.m`
Main script for policy optimization analysis.

**Sections**:
- **MANDATORY PART**: Full dataset (1995-2005) SOP optimization
  - Data import and model setup
  - SOP simulation and visualization
  - NSGA-II optimization
  - Interesting solutions extraction and comparison

- **NON MANDATORY PART**: Calibration/Validation split
  - Dataset splitting (2002 boundary)
  - SOP optimization on calibration
  - Validation performance assessment
  - RBF policy optimization
  - RBF vs SOP comparison

#### `FILE_CONCEGNA_FINALE.m`
Complete ANN forecast model development.

**Sections**:
1. Data loading and deseasonalization
2. K-fold cross-validation partitioning (K=7)
3. Optimal neuron number search (1-5 neurons tested)
4. Optimal lag identification per variable (max_lag=7)
5. Final model training
6. Error analysis (histogram + correlogram)

### Utility Functions

#### `correlogram.m`
```matlab
corr_coeffs = correlogram(X, Y, order, is_unbiased)
```
Computes correlation coefficients between time-lagged signal pairs.

**Parameters**:
- `X`: First signal vector
- `Y`: Second signal vector (default: X for autocorrelation)
- `order`: Maximum lag to compute (default: length(X)-2)
- `is_unbiased`: Boolean for unbiased estimation (default: false)

**Returns**:
- `corr_coeffs`: Vector of correlation coefficients from lag 0 to order

**Features**:
- Automatic plotting if no output argument specified
- 95% confidence bands visualization
- Useful for residual analysis

#### `destagionalizzo.m`
```matlab
[x, ym, sd] = destagionalizzo(u, T, f)
```
Deseasonalizes time series data using moving averages.

**Parameters**:
- `u`: Input time series
- `T`: Period (e.g., 365 for daily data)
- `f`: Half-width of moving average window

**Returns**:
- `x`: Deseasonalized series (standardized)
- `ym`: Cyclo-stationary mean (repeated over series)
- `sd`: Cyclo-stationary std deviation (repeated over series)

**Algorithm**:
1. Compute moving average of original series → mean pattern
2. Compute moving average of squared residuals → variance pattern
3. Standardize: `x = (u - ym) / sd`

---

## 📈 Future Improvements

### Forecast Model
- **Recurrent Neural Networks**: LSTM/GRU to better capture temporal dependencies and reduce error autocorrelation
- **Piecewise Linear Models**: Separate models for dry/wet seasons
- **Ensemble Methods**: Combine multiple model predictions
- **Extended Inputs**: Additional meteorological variables (temperature, humidity)
- **Multi-step Forecasting**: Explicit 5-day ahead prediction vs cumulative

### Policy Optimization
- **Exogenous Information**: Incorporate forecast model into policy (feedforward control)
- **Multi-reservoir System**: Extend to other dams in the basin
- **Dynamic Policies**: Time-varying policies accounting for season
- **Robust Optimization**: Account for forecast uncertainty
- **Additional Objectives**: Water supply, navigation, environmental flows
- **Online Control**: OLFC or POLFC with real-time forecast updates

### Validation
- **Longer Time Series**: More robust solutions with extended data
- **Climate Change Scenarios**: Test policy performance under future conditions
- **Extreme Events**: Specific analysis of rare flood/drought events
- **Sensitivity Analysis**: Parameter robustness assessment

---

## 📚 Theoretical Background

### Standard Operating Policy (SOP)
Piecewise linear control law: `release = f(level)`

**Parameters** (5):
- `h1`: Minimum level to start releasing
- `h2`: Maximum level before spillway activation
- `m1`, `m2`: Slopes of piecewise segments
- `w`: Constant water demand

**Advantages**: Simple, interpretable, commonly used
**Limitations**: Linear, state-only feedback

### Radial Basis Functions (RBF)
Nonlinear control law as weighted sum of Gaussian basis functions:

```
release = scale_factor × Σᵢ wᵢ × φᵢ(level)
φᵢ(level) = exp(-((level - μᵢ)/σᵢ)²)
```

**Parameters** (13):
- 4 means (μ): centers of Gaussian bases
- 4 std deviations (σ): widths of bases
- 4 weights (w): contribution of each base
- 1 scale factor: output scaling

**Advantages**: Flexible, smooth, multivariable-ready
**Limitations**: More parameters, harder to interpret

### NSGA-II Algorithm
Evolutionary multi-objective optimization algorithm.

**Key Features**:
- Fast non-dominated sorting
- Crowding distance for diversity
- Tournament selection
- Polynomial mutation and crossover

**Settings Used**:
- Population: 70 individuals
- Generations: 50
- Decision variables: 5 (SOP) or 13 (RBF)
- Objectives: 2 (hydropower, flooding)

### K-Fold Cross-Validation
Prevents overfitting in time series forecasting.

**Configuration**:
- K = 7 folds (total 7 years)
- Each fold: 5 years calibration, 2 years validation
- Validation years always contiguous
- Progressive shifting of calibration/validation windows

**Metrics**:
- R²: Coefficient of determination
- Mean and std of R² across folds
- Selected configuration: max mean R² in validation

---

## 👥 Authors

**Group 26**:
- Gabriele Ferrari (996460)
- Tommaso Zaini (970230)
- Daniele Sala (996440)

**Course**: Natural Resources Management (A.Y. 2021/2022)  
**Institution**: Politecnico di Milano  
**Program**: Master of Environmental and Land Planning Engineering

**Professors**:
- Andrea Castelletti
- Matteo Giuliani
- Matteo Sangiorgio
- Michele Mauri

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- Politecnico di Milano for providing the case study framework
- Course professors for guidance and support
- Vietnam Electricity (EVN) for the Hoa Binh dam data
- Red River Basin research community

---

## 📞 Contact

For questions or suggestions, please open an issue in this repository or contact the authors.

---

## 📖 References

1. Castelletti, A., Galelli, S., Restelli, M., & Soncini-Sessa, R. (2010). Tree-based reinforcement learning for optimal water reservoir operation. *Water Resources Research*, 46(9).

2. Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. A. M. T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on Evolutionary Computation*, 6(2), 182-197.

3. Giuliani, M., Castelletti, A., Pianosi, F., Mason, E., & Reed, P. M. (2016). Curses, tradeoffs, and scalable management: Advancing evolutionary multiobjective direct policy search to improve water reservoir operations. *Journal of Water Resources Planning and Management*, 142(2).

---

## 🏷️ Keywords

`water-resources-management` `reservoir-optimization` `hydropower` `flood-control` `multi-objective-optimization` `NSGA-II` `EMODPS` `machine-learning` `ANN` `MATLAB` `vietnam` `red-river-basin` `dam-management` `time-series-forecasting` `radial-basis-functions`

---

**Last Updated**: November 2024
