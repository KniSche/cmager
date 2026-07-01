import subprocess
import sys

def verify_r_dependencies():
    """Ensures required R libraries including HGNChelper are present."""
    required_packages = ["readr", "gamsel", "Matrix", "matrixStats", "HGNChelper"]
    
    for pkg in required_packages:
        # Check if the package can be loaded in R
        check_cmd = f"suppressPackageStartupMessages(library({pkg}))"
        result = subprocess.run(["Rscript", "-e", check_cmd], capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"📦 Missing required R package: {pkg}. Attempting auto-installation...")
            
            # The type="binary" directive forces R to use pre-built packages if available on Win/Mac
            install_cmd = (
                f"install.packages('{pkg}', "
                f"repos='https://cloud.r-project.org', "
                f"type='binary')"
            )
            
            try:
                subprocess.run(["Rscript", "-e", install_cmd], check=True)
                print(f"✅ Successfully installed {pkg}!")
            except subprocess.CalledProcessError:
                # If binary isn't an option and it fails on Linux, give a clear platform hint:
                print(
                    f"\n❌ Failed to compile {pkg}.\n"
                    f"💡 This usually means your OS is missing system header files.\n"
                    f"👉 Fix on Ubuntu/Debian: sudo apt install zlib1g-dev libfontconfig1-dev\n"
                    f"👉 Fix on Fedora/CentOS: sudo dnf install zlib-devel\n"
                    f"👉 Fix on macOS: brew install zlib\n",
                    file=sys.stderr
                )
                sys.exit(1)
