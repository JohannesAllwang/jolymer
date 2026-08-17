$ErrorActionPreference = "Stop"

Write-Host "=== Building Jolymer for Windows ==="
Write-Host ""

# ------------------------------------------------------------
# Python environment
# ------------------------------------------------------------

if (Get-Command pyenv -ErrorAction SilentlyContinue) {

    Write-Host "Using pyenv-win"

    $versions = pyenv versions --bare

    if ($versions -notcontains "3.11.9") {
        Write-Host "Installing Python 3.11.9..."
        pyenv install 3.11.9
    }

    if ($versions -notcontains "jolymer") {
        Write-Host "Creating virtual environment 'jolymer'..."
        pyenv virtualenv 3.11.9 jolymer
    }

    pyenv activate jolymer
}

if (-not (Get-Command python -ErrorAction SilentlyContinue)) {
    throw "Python 3.11+ is required."
}

Write-Host "Python version:"
python --version

# ------------------------------------------------------------
# Install package and build dependencies
# ------------------------------------------------------------

Write-Host ""
Write-Host "Installing dependencies..."

python -m pip install --upgrade pip
python -m pip install -e .
python -m pip install nuitka

# ------------------------------------------------------------
# Build
# ------------------------------------------------------------

Write-Host ""
Write-Host "Building standalone application..."

if (Test-Path "build") {
    Remove-Item -Recurse -Force "build"
}

python -m nuitka `
    --standalone `
    --enable-plugin=pyqt6 `
    --include-package=jolymer `
    --output-dir=build `
    --output-filename=jolymer-secplot.exe `
    src/jolymer/gui/sec_plot/app.py

# ------------------------------------------------------------
# Install application
# ------------------------------------------------------------

$installDir = "$env:LOCALAPPDATA\Jolymer"
$distDir = "build\app.dist"

if (Test-Path $installDir) {
    Remove-Item -Recurse -Force $installDir
}

New-Item -ItemType Directory -Force -Path $installDir | Out-Null

Copy-Item `
    -Path "$distDir\*" `
    -Destination $installDir `
    -Recurse

$exe = "$installDir\jolymer-secplot.exe"

# ------------------------------------------------------------
# Create shortcut
# ------------------------------------------------------------

Write-Host ""
$workingDir = Read-Host "Enter the default Jolymer data directory"

if (-not (Test-Path $workingDir)) {
    Write-Host "Directory does not exist. Creating it..."
    New-Item -ItemType Directory -Force -Path $workingDir | Out-Null
}

$shortcutPath = "$env:USERPROFILE\Desktop\Jolymer SEC Plot.lnk"

$WshShell = New-Object -ComObject WScript.Shell
$Shortcut = $WshShell.CreateShortcut($shortcutPath)

$Shortcut.TargetPath = $exe
$Shortcut.WorkingDirectory = $workingDir
$Shortcut.Description = "Jolymer SEC-SWAXS analysis"

$Shortcut.Save()

# ------------------------------------------------------------
# Done
# ------------------------------------------------------------

Write-Host ""
Write-Host "=== Build complete ==="
Write-Host ""
Write-Host "Application:"
Write-Host "    $exe"
Write-Host ""
Write-Host "Shortcut:"
Write-Host "    $shortcutPath"
Write-Host ""
