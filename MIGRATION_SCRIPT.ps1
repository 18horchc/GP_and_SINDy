# Safe Migration Script for GP_and_SINDy Project
# This script moves the git repository from nested structure to clean structure
# 
# BEFORE RUNNING:
# 1. Ensure all changes are committed
# 2. Close MATLAB and all editors
# 3. Create a backup of the entire folder
# 4. Review this script carefully

Write-Host "=== GP_and_SINDy Structure Migration ===" -ForegroundColor Cyan
Write-Host ""

# Step 1: Verify we're in the right location
$currentDir = Get-Location
Write-Host "Current directory: $currentDir" -ForegroundColor Yellow

if (-not (Test-Path "GP_and_SINDy\.git")) {
    Write-Host "ERROR: Git repository not found in expected location!" -ForegroundColor Red
    Write-Host "Expected: $currentDir\GP_and_SINDy\.git" -ForegroundColor Red
    exit 1
}

Write-Host "✓ Git repository found" -ForegroundColor Green

# Step 2: Verify git status is clean
Write-Host ""
Write-Host "Checking git status..." -ForegroundColor Yellow
Push-Location "GP_and_SINDy"
$gitStatus = git status --porcelain
if ($gitStatus) {
    Write-Host "WARNING: Uncommitted changes detected!" -ForegroundColor Red
    Write-Host $gitStatus
    Write-Host ""
    $response = Read-Host "Continue anyway? (y/n)"
    if ($response -ne "y") {
        Pop-Location
        exit 1
    }
} else {
    Write-Host "✓ Git status is clean" -ForegroundColor Green
}
Pop-Location

# Step 3: Check for processes using files
Write-Host ""
Write-Host "IMPORTANT: Please ensure:" -ForegroundColor Yellow
Write-Host "  1. MATLAB is closed" -ForegroundColor Yellow
Write-Host "  2. All editors are closed" -ForegroundColor Yellow
Write-Host "  3. No file explorers are open in this directory" -ForegroundColor Yellow
Write-Host ""
$response = Read-Host "Are all applications closed? (y/n)"
if ($response -ne "y") {
    Write-Host "Please close all applications and run this script again." -ForegroundColor Red
    exit 1
}

# Step 4: Create backup reminder
Write-Host ""
Write-Host "BACKUP REMINDER:" -ForegroundColor Yellow
Write-Host "Have you created a backup? (Recommended: copy entire GP_and_SINDy folder)" -ForegroundColor Yellow
$response = Read-Host "Continue with migration? (y/n)"
if ($response -ne "y") {
    Write-Host "Migration cancelled." -ForegroundColor Yellow
    exit 0
}

# Step 5: Perform migration
Write-Host ""
Write-Host "Starting migration..." -ForegroundColor Cyan

try {
    # Move .git folder
    Write-Host "Moving .git folder..." -ForegroundColor Yellow
    Move-Item -Path "GP_and_SINDy\.git" -Destination ".git" -Force
    Write-Host "✓ .git moved" -ForegroundColor Green
    
    # Move .gitignore
    if (Test-Path "GP_and_SINDy\.gitignore") {
        Write-Host "Moving .gitignore..." -ForegroundColor Yellow
        Move-Item -Path "GP_and_SINDy\.gitignore" -Destination ".gitignore" -Force
        Write-Host "✓ .gitignore moved" -ForegroundColor Green
    }
    
    # Move .gitattributes if exists
    if (Test-Path "GP_and_SINDy\.gitattributes") {
        Write-Host "Moving .gitattributes..." -ForegroundColor Yellow
        Move-Item -Path "GP_and_SINDy\.gitattributes" -Destination ".gitattributes" -Force
        Write-Host "✓ .gitattributes moved" -ForegroundColor Green
    }
    
    # Move all files and folders from inner to outer
    Write-Host "Moving all files and folders..." -ForegroundColor Yellow
    Get-ChildItem -Path "GP_and_SINDy" -Force | ForEach-Object {
        $dest = Join-Path "." $_.Name
        if (Test-Path $dest) {
            Write-Host "  Skipping $($_.Name) - already exists in destination" -ForegroundColor Yellow
        } else {
            Move-Item -Path $_.FullName -Destination $dest -Force
            Write-Host "  Moved $($_.Name)" -ForegroundColor Gray
        }
    }
    Write-Host "✓ All files moved" -ForegroundColor Green
    
    # Remove empty inner folder
    Write-Host "Removing empty inner folder..." -ForegroundColor Yellow
    Remove-Item -Path "GP_and_SINDy" -Force -ErrorAction SilentlyContinue
    Write-Host "✓ Inner folder removed" -ForegroundColor Green
    
} catch {
    Write-Host "ERROR during migration: $_" -ForegroundColor Red
    Write-Host "Please check the state of your repository!" -ForegroundColor Red
    exit 1
}

# Step 6: Verify migration
Write-Host ""
Write-Host "Verifying migration..." -ForegroundColor Cyan
Push-Location "."

# Check git still works
Write-Host "Testing git..." -ForegroundColor Yellow
$gitTest = git status
if ($LASTEXITCODE -eq 0) {
    Write-Host "✓ Git is working" -ForegroundColor Green
} else {
    Write-Host "ERROR: Git is not working!" -ForegroundColor Red
    Pop-Location
    exit 1
}

# Check git log
Write-Host "Verifying git history..." -ForegroundColor Yellow
$gitLog = git log --oneline -1
if ($gitLog) {
    Write-Host "✓ Git history intact: $gitLog" -ForegroundColor Green
} else {
    Write-Host "WARNING: Could not verify git history" -ForegroundColor Yellow
}

# Check remote
Write-Host "Verifying remote..." -ForegroundColor Yellow
$gitRemote = git remote -v
if ($gitRemote) {
    Write-Host "✓ Remote configured:" -ForegroundColor Green
    Write-Host $gitRemote
} else {
    Write-Host "WARNING: No remote configured" -ForegroundColor Yellow
}

Pop-Location

# Step 7: Final instructions
Write-Host ""
Write-Host "=== Migration Complete ===" -ForegroundColor Green
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "1. Test git operations: git status, git log" -ForegroundColor White
Write-Host "2. Test remote: git fetch (don't push yet)" -ForegroundColor White
Write-Host "3. Open MATLAB project: GP_and_SINDy.prj" -ForegroundColor White
Write-Host "4. Verify all files are accessible in MATLAB" -ForegroundColor White
Write-Host "5. Run a test script to verify paths work" -ForegroundColor White
Write-Host "6. If everything works, you can push: git push" -ForegroundColor White
Write-Host ""
Write-Host "Current directory is now the git root: $(Get-Location)" -ForegroundColor Green

