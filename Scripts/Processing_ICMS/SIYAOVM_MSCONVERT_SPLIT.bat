@echo off
setlocal enabledelayedexpansion

REM --------------------------------------------
REM >>> EDIT THIS LINE IF NEEDED <<<
set MSCONVERT_EXE=C:\Program Files\ProteoWizard\ProteoWizard 3.0.24066.8d99481\msconvert.exe

REM If MSCONVERT_EXE not set, assume it's in PATH
if "%MSCONVERT_EXE%"=="" (
    set MSCONVERT_EXE=msconvert.exe
)

echo Using msconvert: "%MSCONVERT_EXE%"

REM Create output folder "LM-mzml" if not exists
if not exist "LM-mzml" mkdir "LM-mzml"

REM Loop through all .raw files in current folder
for %%F in (*.raw) do (
    echo Processing LM version of %%F ...

    REM Get creation date in YYYYMMDD format
    for /f "usebackq tokens=*" %%D in (`powershell -NoProfile -Command ^
        "(Get-Item '%%~fF').CreationTime.ToString('yyyyMMdd')"`) do set "filedate=%%D"

    REM Run msconvert with LM filters and date+prefix in output filename
    call "%MSCONVERT_EXE%" "%%F" ^
        --zlib ^
        --filter "peakPicking vendor msLevel=1-" ^
        --filter "thermoScanFilter contains include 65.0000-250.0000" ^
        --outdir "LM-mzml" ^
        --outfile "!filedate!_LM-zHILIC-IQX_LmexQC_%%~nF.mzML"

    REM Copy timestamps from raw to mzML
    set "mzmlfile=LM-mzml\!filedate!_LM-zHILIC-IQX_LmexQC_%%~nF.mzML"
    powershell -Command ^
      "$src=Get-Item '%%~fF'; $dst=Get-Item '!mzmlfile!';" ^
      "$dst.CreationTime=$src.CreationTime; $dst.LastWriteTime=$src.LastWriteTime"
)

REM Create output folder "HM-mzml" if not exists
if not exist "HM-mzml" mkdir "HM-mzml"

REM Loop through all .raw files in current folder
for %%F in (*.raw) do (
    echo Processing HM version of %%F ...

    REM Get creation date in YYYYMMDD format
    for /f "usebackq tokens=*" %%D in (`powershell -NoProfile -Command ^
        "(Get-Item '%%~fF').CreationTime.ToString('yyyyMMdd')"`) do set "filedate=%%D"

    REM Run msconvert with HM filters and date+prefix in output filename
    call "%MSCONVERT_EXE%" "%%F" ^
        --zlib ^
        --filter "peakPicking vendor msLevel=1-" ^
        --filter "thermoScanFilter contains include 195.0000-1900.0000" ^
        --outdir "HM-mzml" ^
        --outfile "!filedate!_HM-zHILIC-IQX_LmexQC_%%~nF.mzML"

    REM Copy timestamps from raw to mzML
    set "mzmlfile=HM-mzml\!filedate!_HM-zHILIC-IQX_LmexQC_%%~nF.mzML"
    powershell -Command ^
      "$src=Get-Item '%%~fF'; $dst=Get-Item '!mzmlfile!';" ^
      "$dst.CreationTime=$src.CreationTime; $dst.LastWriteTime=$src.LastWriteTime"
)

echo All conversions finished!
pause
