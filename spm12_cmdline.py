import subprocess
import sys
if sys.version_info[0] < 3:
    raise Exception("Python 3.0+ is needed.")
    
subprocess.run(["python", "./spm12_scripts/spm12_cmdline_main.py"] + sys.argv[1:])
