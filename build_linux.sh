python3 -m nuitka \
--standalone \
--enable-plugin=numpy \
--include-data-dir=sturgeon/include=sturgeon/include \
sturgeon/main.py 

python3 -m nuitka \
--onefile \
--enable-plugin=numpy \
--include-data-dir=sturgeon/include=sturgeon/include \
sturgeon/main.py 
