Repository for refactoring the "MAWS 2023" repository by `dtu-denmark` for integration in the [pyaptamer](https://github.com/gc-os-ai/pyaptamer/) package

### Refactoring aims

* remove or replace binary dependencies
* remove or replace subprocess calls
* full python bindings

### Local development & testing 🛠️

To reproduce the environment used in this workspace, run the following commands from a
shell (Git Bash, WSL, PowerShell with appropriate tools, etc.):

```bash
# 1. clone the repository (if not already done)
git clone https://github.com/gc-os-ai/MAWS_2025.git
cd MAWS_2025

# 2. (preferred) create a conda environment from the provided YAML file
conda env create -f environment.yml
conda activate maws

# alternatively, you can use a plain venv:
# python -m venv .venv
# Windows PowerShell:
# . .venv/Scripts/Activate.ps1
# bash / WSL:
# source .venv/bin/activate
# pip install -r requirements.txt
# pip install openmm numba

# 3. install any additional tools
pip install pytest nbformat nbconvert jupyter pre-commit

# 4. run the test suite (includes notebook execution)
python -m pytest -q

# 5. run pre-commit manually if you want to check formatting
pre-commit run --all-files
```

The tests added in `tests/` verify that all modules import cleanly and that basic
utility functions behave correctly. The examples in `examples/` can be executed as a
higher-level sanity check but they depend on AmberTools (`tleap`, `antechamber`,
`parmchk2`) being available on your `PATH`.


### original readme

https://github.com/gc-os-ai/MAWS_2025/blob/main/README_orig.md
