# python3 -m nuitka \
# --standalone \
# --enable-plugin=numpy \
# --include-data-dir=sturgeon/include=sturgeon/include \
# sturgeon/main.py 

python3 -m nuitka \
--onefile \
--enable-plugin=numpy \
--include-data-dir=sturgeon/include=sturgeon/include \
sturgeon/main.py 

rm bin/sturgeon_macos12_86_64.bin
mv main.bin bin/sturgeon_macos12_86_64.bin