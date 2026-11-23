export NTE_ROOT=/home/ps/jishaobin/research-field-NTE

export CONDARC="$NTE_ROOT/condarc-local.yml"
export CONDA_PKGS_DIRS="$NTE_ROOT/.cache/conda-pkgs"
export CONDA_ENVS_PATH="$NTE_ROOT/envs"
export CONDA_REGISTER_ENVS=false

export XDG_CACHE_HOME="$NTE_ROOT/.cache"
export PIP_CACHE_DIR="$NTE_ROOT/.cache/pip"
export MPLCONFIGDIR="$NTE_ROOT/.cache/matplotlib"
export CUDA_CACHE_PATH="$NTE_ROOT/.cache/cuda"
export TORCH_HOME="$NTE_ROOT/.cache/torch"
export TMPDIR="$NTE_ROOT/.tmp"

## https://github.com/vllm-project/vllm/issues/10300
export LD_LIBRARY_PATH=$NTE_ROOT/envs/deepsdf/lib/python3.10/site-packages/nvidia/nvjitlink/lib:$LD_LIBRARY_PATH

export PYTHONNOUSERSITE=1
export PIP_CONFIG_FILE=/dev/null
export PIP_USER=0
export MPLBACKEND=Agg

unset PYTHONPATH PYTHONHOME PIP_TARGET PIP_PREFIX
