
cmake . \
    -DEDRIXS_PY_INTERFACE=ON \
    -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=build/lib.linux-x86_64-cpython-39/edrixs\
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=build/temp.linux-x86_64-cpython-39 \
    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=build/temp.linux-x86_64-cpython-39 \
    -DCMAKE_LIBRARY_PATH=$HOME/.local/lib \
