@echo off
cd /d "%~dp0"
python run_table4.py
if errorlevel 1 (
    echo.
    echo Table 4 failed.
    pause
    exit /b 1
)
pause
