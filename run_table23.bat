@echo off
cd /d "%~dp0"
python run_table23.py
if errorlevel 1 (
    echo.
    echo Tables 2-3 failed.
    pause
    exit /b 1
)
pause
