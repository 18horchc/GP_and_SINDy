# Migration Checklist - Action Required

## ✅ Completed
- [x] All changes committed to git
- [x] Migration script created
- [x] Risk analysis documented

## ⚠️ Pre-Migration Actions Required (YOUR ACTION NEEDED)

### 1. Close All Applications
**Please close:**
- [ ] MATLAB (completely closed, not just minimized)
- [ ] VS Code / Cursor / any code editors
- [ ] File Explorer windows in the GP_and_SINDy folder
- [ ] Any other programs that might have files open

**How to verify:** Check Task Manager - no MATLAB.exe or editor processes running

### 2. Create Backup (STRONGLY RECOMMENDED)
**Please create a backup:**
- [ ] Copy the entire `C:\Users\chorc\GP_and_SINDy` folder to a safe location
- [ ] Suggested location: `C:\Users\chorc\GP_and_SINDy_BACKUP_YYYYMMDD`
- [ ] Verify the backup was created successfully

### 3. Duplicate Files Found
I found these files in BOTH outer and inner folders. They appear to be duplicates, but I need your confirmation:

**Files in outer folder that also exist in inner:**
- `build_library.m`
- `degrade_data.m`
- `fitAndPlotGP.m`
- `lotka_volterra_comparison.m`
- `LOTKA_VOLTERRA_README.md`
- `Lotka_Volterra_toy_problem.m`
- `LV_ground_truth.m`
- `LV_toy_problem_2.m`
- `README.md`
- `run_ESINDy.m`
- `runGP.m`

**Question:** Are these files in the outer folder:
- A) Old duplicates that can be safely deleted?
- B) Different versions that need to be compared/merged?
- C) You're not sure - need me to compare them?

## 📋 Migration Plan

Once you've completed the checklist above, I will:

1. **Compare duplicate files** (if needed) to ensure no data loss
2. **Run the migration script** which will:
   - Move `.git` folder to outer directory
   - Move all files from inner to outer
   - Remove empty inner folder
   - Verify git still works
3. **Test thoroughly**:
   - Verify git status, log, remote
   - Check file integrity
   - Test MATLAB project access

## 🎯 Ready to Proceed?

Once you've:
- ✅ Closed all applications
- ✅ Created a backup
- ✅ Answered the duplicate files question

**Reply with:** "Ready to migrate" and I'll proceed with the safe migration.

