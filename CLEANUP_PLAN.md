# Project Structure Cleanup Plan

## Current Situation
- **Git Repository**: Located in `C:\Users\chorc\GP_and_SINDy\GP_and_SINDy\` (inner folder)
- **Outer Folder**: `C:\Users\chorc\GP_and_SINDy\` contains some duplicate files
- **Issue**: Nested structure causing MATLAB project visibility problems

## Unique Files in Outer Folder
- `gp_microglia_fig11_demo.m` (needs to be preserved)
- `gp_microglia_fig11_demo.asv` (autosave, can be ignored)

## Recommended Solution: Move Git Repo to Outer Folder

**Why this approach:**
- Cleaner structure (no nested folders)
- MATLAB project can easily reference all files
- Matches typical project structure
- Easier to maintain

**Steps:**
1. Copy unique file from outer to inner (backup)
2. Move .git folder from inner to outer
3. Move all files from inner to outer
4. Update git config if needed
5. Remove empty inner folder
6. Update MATLAB project file

## Alternative: Keep Git in Inner, Update MATLAB Project

**Why this approach:**
- Safer (less moving of git)
- Minimal changes

**Steps:**
1. Copy unique file from outer to inner
2. Update MATLAB project file to point to inner folder
3. Remove duplicate files from outer folder
4. Document structure

## Recommendation
I recommend **Option 1** (move git to outer) for cleaner structure, but we'll do it safely with backups.

