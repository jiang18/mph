## Source code
The `mph` source code is hosted on [Github](https://github.com/jiang18/mph).

## Installation
### Binary executables
- Binary executables are available in the [Releases](https://github.com/jiang18/mph/releases/latest).
- `mph` is statically compiled with [Intel oneAPI Base & HPC Toolkit 2023.1](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html) and Eigen 3.4.0 on Ubuntu 22.
- `mph` should be executable. Type `chmod +x mph` on the command line if needed.

### Compilation of `mph`
1. Requirements
    - Linux x86_64
    - [Intel oneAPI Base & HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html)
    - Eigen >= 3.4.0
2. Change the Eigen path in `Makefile`. 
3. On the command line type `make` while in the directory containing the `Makefile`.
4. This should produce the executable named `mph`.

## For Windows
Windows Subsystem for Linux (WSL) can be used to run MPH on a Windows machine.

1. [Install Linux on Windows with WSL](https://learn.microsoft.com/en-us/windows/wsl/install)
    - Ubuntu will be installed by default, and other Linux distributions are available. 
2. Run [pre-compiled MPH binaries](https://github.com/jiang18/mph/releases/latest) on a WSL Linux distribution as on a Linux machine.