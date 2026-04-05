param(
  [int]$L = 32,
  [int]$NTherm = 3000,
  [int]$NSweep = 20,
  [int]$NMc = 5000,
  [string]$PythonExe = "python"
)

$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

Write-Host "[1/3] Build ising.x ..."
if (Get-Command g++ -ErrorAction SilentlyContinue) {
  g++ -O3 -o ising.x main.C ising.C mc.C measurement.C mt19937ar_lib.C stat_utils.C
} elseif (Get-Command clang++ -ErrorAction SilentlyContinue) {
  clang++ -O3 -o ising.x main.C ising.C mc.C measurement.C mt19937ar_lib.C stat_utils.C
} else {
  throw "g++ または clang++ が見つかりません。MSYS2/MinGW などをインストールしてください。"
}

Write-Host "[2/3] Install plotting dependencies ..."
& $PythonExe -m pip install --user matplotlib numpy

Write-Host "[3/3] Run simulation + plot ..."
& $PythonExe .\plot_phase_transition.py --exe .\ising.x --L $L --n-therm $NTherm --n-sweep $NSweep --n-mc $NMc

Write-Host "Done. Check .\results\*.png"
