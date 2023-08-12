_run FILE:
    ./{{FILE}}
_debug FILE:
    valgrind ./{{FILE}}
@run-fast FILE:
    cc {{FILE}}.c rlFPCamera/rlFPCamera.c -O3 -Wall -Wextra -ffast-math -funsafe-math-optimizations -fno-math-errno -l:libraylib.a -lm -pthread -o {{FILE}}
    just _run {{FILE}}
@compile FILE:
    cc {{FILE}}.c rlFPCamera/rlFPCamera.c -Og -Wall -Wextra -fanalyzer -l:libraylib.a -lm -pthread -o {{FILE}}
@run FILE:
    cc {{FILE}}.c rlFPCamera/rlFPCamera.c -Og -Wall -Wextra -l:libraylib.a -lm -pthread -o {{FILE}}
    just _run {{FILE}}
@run_web FILE:
    emcc {{FILE}}.c rlFPCamera/rlFPCamera.c -Os -Wall -Wextra -DPLATFORM_WEB -sASSERTIONS -sALLOW_MEMORY_GROWTH -I/usr/local/include -L/usr/local/lib ./libraylib.a --emrun -s USE_GLFW=3 --shell-file shell_minimal.html -o {{FILE}}.html 
    emrun {{FILE}}.html
@compile_web FILE:
    emcc {{FILE}}.c rlFPCamera/rlFPCamera.c -Os -Wall -Wextra -DPLATFORM_WEB -sASSERTIONS -sALLOW_MEMORY_GROWTH -I/usr/local/include -L/usr/local/lib ./libraylib.a  -s USE_GLFW=3 --shell-file shell_minimal.html -o {{FILE}}.html 
@debug FILE:
    cc {{FILE}}.c rlFPCamera/rlFPCamera.c -Og -g3 -Wall -Wextra -l:libraylib.a -lm -pthread -o {{FILE}}
    just _debug {{FILE}}
@format FILE:
    clang-format -style=file:.llvm-format -i {{FILE}}.c 
