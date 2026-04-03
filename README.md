# DA2MODE-Algorithm
DA2MODE: Dynamic Archive with Adaptive Multi-Operator Differential Evolution for numerical optimization

Official MATLAB implementation of the DA2MODE algorithm.

## 📄 Paper

Mohamed Reda, Ahmed Onsy, Amira Y. Haikal, and Ali Ghanbari  
**DA2MODE: Dynamic Archive with Adaptive Multi-Operator Differential Evolution for numerical optimization**  
Swarm and Evolutionary Computation, Volume 99, 2025, Article 102130  
DOI: https://doi.org/10.1016/j.swevo.2025.102130

---

## Overview

DA2MODE is a Differential Evolution (DE) algorithm designed to improve:

- Adaptive operator selection
- Population diversity
- Convergence speed
- Archive-driven exploration

Main contributions:
- Progressive Adaptive Selector with Exponential Smoothing (PASES)
- Adaptive Non-Elite Archive Update (ANEAU)
- Adaptive mutation and crossover selection
- Adaptive control parameters

---

## 📂 Repository Structure

```text
DA2MODE-Algorithm/
│
├── DA2MODE_algorithm.m
├── run_DA2MODE_main.m
├── cost_cec2017.m
├── cost_cec2020.m
├── cost_cec2022.m
│
├── cec17_func.cpp
├── cec17_func.mex*
├── cec20_func.cpp
├── cec20_func.mex*
├── cec22_test_func.cpp
├── cec22_test_func.mex*
│
├── input_data2017/
├── input_data2020/
├── input_data2022/
│
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
└── .gitattributes
```

---

## ⚙️ Requirements

- MATLAB (recommended R2023a or later)
- Statistics and Machine Learning Toolbox (`lhsdesign`, `normrnd`)

---

## ▶️ How to Run

You can run the algorithm in three ways.

### Option 1: Run with default values
If no input arguments are passed, the algorithm uses its internal default settings.

```matlab
[goalReached, GlobalBest, countFE] = DA2MODE_algorithm();
```

Default values inside `DA2MODE_algorithm.m`:

```matlab
CECyear = 2020;
fNo     = 3;
nd      = 20;
lb      = -100;
ub      = 100;
```

### Option 2: Run by passing parameters directly

```matlab
[goalReached, GlobalBest, countFE] = DA2MODE_algorithm(CECyear, fNo, nd, lb, ub);
```

Example for CEC2020:

```matlab
[goalReached, GlobalBest, countFE] = DA2MODE_algorithm(2020, 3, 20, -100, 100);
```

Example for CEC2022:

```matlab
[goalReached, GlobalBest, countFE] = DA2MODE_algorithm(2022, 5, 10, -100, 100);
```

Example for CEC2017:

```matlab
[goalReached, GlobalBest, countFE] = DA2MODE_algorithm(2017, 10, 30, -100, 100);
```

### Option 3: Run using the main configuration file

```matlab
run_DA2MODE_main
```

---

## 🔧 Input Parameters

The function format is:

```matlab
DA2MODE_algorithm(CECyear, fNo, nd, lb, ub)
```

### Parameters Description

- `CECyear` : Benchmark year  
  Supported values:
  - `2017`
  - `2020`
  - `2022`

---

- `fNo` : Benchmark function number  

  - For **CEC2017**:
    - valid range: `1 – 29`

  - For **CEC2020**:
    - valid range: `1 – 10`

  - For **CEC2022**:
    - valid range: `1 – 12`

---

- `nd` : Problem dimension  

  Supported values in this implementation:
  - For **CEC2020**: `10`, `20`
  - For **CEC2022**: `10`, `20`
  - For **CEC2017**: `30`, `50`, `100`

---

- `lb` : Lower bound  
  Typically:
  ```matlab
  -100
  ```

- `ub` : Upper bound  
  Typically:
  ```matlab
  100
  ```

---

## 📌 Important Notes

- Make sure that:
  - the correct CEC input data folders are available
  - the corresponding benchmark functions are compiled as MEX files
  - the benchmark year, function number, and dimension match the selected CEC suite

- The benchmark search space is typically:
  ```matlab
  [-100, 100]
  ```

- The stopping criterion in the published paper uses a tolerance of:
  ```matlab
  1e-8
  ```

---

## 🧪 Example Runs

### CEC2020 (Function 3, 20D)
```matlab
DA2MODE_algorithm(2020, 3, 20, -100, 100);
```

### CEC2020 (Function 10, 10D)
```matlab
DA2MODE_algorithm(2020, 10, 10, -100, 100);
```

### CEC2022 (Function 5, 20D)
```matlab
DA2MODE_algorithm(2022, 5, 20, -100, 100);
```

### CEC2022 (Function 12, 10D)
```matlab
DA2MODE_algorithm(2022, 12, 10, -100, 100);
```

### CEC2017 (Function 8, 30D)
```matlab
DA2MODE_algorithm(2017, 8, 30, -100, 100);
```

### CEC2017 (Function 29, 100D)
```matlab
DA2MODE_algorithm(2017, 29, 100, -100, 100);
```

---

## 📚 Citation

If you use this code, please cite:

```bibtex
@article{reda2025da2mode,
  title   = {DA2MODE: Dynamic Archive with Adaptive Multi-Operator Differential Evolution for numerical optimization},
  author  = {Reda, Mohamed and Onsy, Ahmed and Haikal, Amira Y. and Ghanbari, Ali},
  journal = {Swarm and Evolutionary Computation},
  volume  = {99},
  pages   = {102130},
  year    = {2025},
  publisher = {Elsevier},
  doi     = {10.1016/j.swevo.2025.102130}
}
```

---

## 📜 License

This project is released under the MIT License. See the `LICENSE` file for details.

---

## 📧 Contact

**Dr. Mohamed Reda**  
University of Central Lancashire, UK  
Mansoura University, Egypt

- 📩 Personal: [mohamed.reda.mu@gmail.com](mailto:mohamed.reda.mu@gmail.com)  
- 📩 Academic: [mohamed.reda@mans.edu.eg](mailto:mohamed.reda@mans.edu.eg)

---

## 🌐 Academic Profiles

- 🧑‍🔬 ORCID: https://orcid.org/0000-0002-6865-1315  
- 🎓 Google Scholar: https://scholar.google.com/citations?user=JmuB2qwAAAAJ  
- 📊 Scopus: https://www.scopus.com/authid/detail.uri?authorId=57220204540  
- 📚 Web of Science: https://www.webofscience.com/wos/author/record/3164983  
- 🧾 SciProfiles: https://sciprofiles.com/profile/Mreda  

---

## 🔗 Professional & Social Links

- 💼 LinkedIn: https://www.linkedin.com/in/mraf  
- 🔬 ResearchGate: https://www.researchgate.net/profile/Mohamed-Reda-8  
- 🎓 Academia: https://mansoura.academia.edu/MohamedRedaAboelfotohMohamed  
- 📘 SciLit: https://www.scilit.net/scholars/12099081  
- 🧮 MATLAB Central: https://uk.mathworks.com/matlabcentral/profile/authors/36082525  
- ▶️ YouTube: https://youtube.com/@mredacs  

---

## Acknowledgement

This repository accompanies the published DA2MODE paper in *Swarm and Evolutionary Computation*.

