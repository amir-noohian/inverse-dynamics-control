# Inverse Dynamics Control for the Barrett WAM (4-DOF)

Implementations and experiments for **gravity compensation**, **dynamic compensation**, and **inverse dynamics control**  on a 4-DOF Barrett **WAM** robotic arm.

---

## ✨ Features

- **Inverse Dynamics Control** for 4-DOF WAM.
- **Dynamic Compensation** and **Gravity Compensation**.
- Multiple build/run targets for **Ares** and **Zeus** setups.
- Python tools to visualize tracking and compensation performance (position/torque).

---

## 📁 Repository Structure

- `4dof_dynamics_compensation.cpp` – Dynamic compensation (+ friction).
- `4dof_dyncomp_ares.cpp`, `4dof_dyncomp_zeus.cpp` – Dynamic comp entry points per setup.
- `4dof_gravcomp_ares.cpp`, `4dof_gravcomp_zeus.cpp` – Gravity compensation apps.
- `4dof_id_control.cpp` – Inverse dynamics control app.
- `include/`, `config/`, `build/`, `.vscode/`, `scripts/`, `data/`, `old_1/` – Headers, configs, logs/plots, helpers, and archived versions.
- `plot_dynamics_compensation.py`, `plot_gravcomp_{ares,zeus}.py`, `plot_id_control.py` – Plotting utilities.
- `CMakeLists.txt` – Build configuration.

---

## 🛠️ Build

### Prerequisites
- **C++17**, **CMake ≥ 3.10**
- Linear algebra library (e.g., **Eigen**)
- Barrett WAM dependencies (e.g., `libbarrett`) if running on hardware
- (Optional) **Python 3.x** with `numpy`, `matplotlib` for plotting

### Compile
```bash
mkdir -p build && cd build
cmake ..
make -j
```

### Run (examples)
```bash
# Dynamic compensation (Ares / Zeus)
./4dof_dyncomp_ares
./4dof_dyncomp_zeus

# Gravity compensation (Ares / Zeus)
./4dof_gravcomp_ares
./4dof_gravcomp_zeus

# Inverse dynamics control
./4dof_id_control
```

> Notes:
> - Some targets expect configuration files under `config/` and will log to `data/`.
> - Hardware runs require a properly configured WAM and the Barrett stack.

---

## 📊 Plotting & Analysis

Use the Python utilities to visualize logs (position/torque/compensation).

```bash
# Dynamic compensation plots
python3 plot_dynamics_compensation.py

# Gravity compensation plots
python3 plot_gravcomp_ares.py
python3 plot_gravcomp_zeus.py

# Inverse dynamics control plots
python3 plot_id_control.py
```

---

## 🧠 Background

Let \(q,\dot q,\ddot q\) be the joint position, velocity, and acceleration.  
The nominal rigid-body dynamics are:

\[
\tau = M(q)\ddot q + C(q,\dot q)\dot q + g(q) + \tau_f(\dot q)
\]

The controllers are implemented in three stages:

1. **Gravity Compensation** – cancels the gravitational term \(g(q)\).  
2. **Inverse Dynamics Control** – applies the full model \(M(q), C(q,\dot q), g(q)\) for motion tracking.  
3. **Inverse Dynamics with Friction** – extends the model with friction \(\tau_f(\dot q)\) (Coulomb and viscous/stiction), improving accuracy in slow and precise motions.


---

## 📝 Citation

If this repo helps your research, please cite it or link to this page.

---

## 📄 License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.
