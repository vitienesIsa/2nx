@echo off
title H5 to NXtomo Converter

echo ============================================
echo   H5 to NXtomo Converter
echo ============================================
echo.

:: <-- replace "" with your proxy address if needed, e.g. "http://proxy.institution.de:8080"
set PROXY="" 


:: Check Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python was not found. Please install Python 3.13 from https://www.python.org
    echo Make sure to tick "Add Python to PATH" during installation.
    pause
    exit /b 1
)

:: Install/verify dependencies
echo Checking dependencies...
if %PROXY%=="" (
    pip install h5py nxtomo silx pint tqdm numpy --quiet
) else (
    pip install h5py nxtomo silx pint tqdm numpy --quiet --proxy %PROXY%
)
if errorlevel 1 (
    echo.
    echo ERROR: Failed to install dependencies.
    echo Try running this file as Administrator, or install manually with:
    echo   pip install h5py nxtomo silx pint tqdm numpy
    pause
    exit /b 1
)
echo Dependencies OK.
echo.

:: Run the converter
echo Starting converter...
echo.
python "%~dp0h5tonx.py"

if errorlevel 1 (
    echo.
    echo ERROR: The converter exited with an error. See message above.
    pause
    exit /b 1
)

exit /b 0
