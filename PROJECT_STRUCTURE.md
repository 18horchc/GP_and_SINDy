# Project Structure Documentation

## Current Structure

```
C:\Users\chorc\GP_and_SINDy\              (Outer folder)
├── GP_and_SINDy\                         (Inner folder - **GIT REPOSITORY IS HERE**)
│   ├── .git\                             (Git repository)
│   ├── Van_der_Pol\                      (New Van der Pol implementation)
│   ├── Lotka_Volterra\                   (Lotka-Volterra implementations)
│   ├── SarasCode\                        (SINDy code)
│   ├── Archive\                          (Old code)
│   ├── GP_and_SINDy.prj                  (MATLAB project file - **USE THIS ONE**)
│   └── [all other project files]
│
└── [Some duplicate/old files - can be ignored]
    ├── GP_and_SINDy.prj                  (Old project file - ignore)
    └── gp_microglia_fig11_demo.m         (Copied to inner folder)
```

## Important Notes

### ✅ Where to Work
- **All work should be done in**: `C:\Users\chorc\GP_and_SINDy\GP_and_SINDy\`
- **This is where the git repository is located**
- **This is where all tracked files are**

### 📁 MATLAB Project File
- **Use this project file**: `C:\Users\chorc\GP_and_SINDy\GP_and_SINDy\GP_and_SINDy.prj`
- Open this file in MATLAB to have access to all project files
- The project file in the outer folder is outdated and can be ignored

### 🔧 Git Operations
- Always run git commands from: `C:\Users\chorc\GP_and_SINDy\GP_and_SINDy\`
- This is the git repository root

### 📝 Why This Structure Exists
The nested structure likely happened during initial project setup. The git repository was initialized in the inner folder, and some files ended up in the outer folder.

### 🧹 Cleanup Status
- ✅ Unique file `gp_microglia_fig11_demo.m` has been copied to git repo
- ✅ All project files are in the git repository
- ✅ MATLAB project file is in the correct location

### 🚀 How to Use

1. **Open MATLAB**
2. **Open Project**: File → Open Project → Navigate to `C:\Users\chorc\GP_and_SINDy\GP_and_SINDy\GP_and_SINDy.prj`
3. **All files will be accessible** in the MATLAB project interface

### ⚠️ Future Cleanup (Optional)

If you want to eliminate the nested structure entirely, you would need to:
1. Move `.git` folder from inner to outer
2. Move all files from inner to outer
3. Update git config
4. Remove inner folder

**However, this is not necessary** - the current structure works fine as long as you open the project from the inner folder.

