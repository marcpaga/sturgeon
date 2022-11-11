export INCLUDE=/usr/include
export C_INCLUDE_PATH=/usr/include
export CPLUS_INCLUDE_PATH=/usr/include
export LIBRARY_PATH=/usr/lib/x86_64-conda-linux-gnu-ld

python3 -m nuitka \
--standalone \
--enable-plugin=numpy \
--include-data-dir=sturgeon/include=include \
sturgeon/main.py 

python3 -m nuitka \
--onefile \
--enable-plugin=numpy \
--include-data-dir=sturgeon/include=include \
sturgeon/main.py 
