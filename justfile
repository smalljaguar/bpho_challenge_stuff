_run FILE:
    ./{{FILE}}
_debug FILE:
    valgrind ./{{FILE}}
@run-fast FILE:
    cc {{FILE}}.c -O3 -Wall -Wextra -ffast-math -funsafe-math-optimizations -l:libraylib.a -lm -pthread -o {{FILE}}
    just _run {{FILE}}
@run FILE:
    cc {{FILE}}.c -O3 -Wall -Wextra -l:libraylib.a -lm -pthread -o {{FILE}}
    just _run {{FILE}}
@run_web FILE:
    emcc {{FILE}}.c -Os -Wall -Wextra -DPLATFORM_WEB -sASSERTIONS -sALLOW_MEMORY_GROWTH -I/usr/local/include -L/usr/local/lib ./libraylib.a --emrun -s USE_GLFW=3 -o {{FILE}}.html 
    emrun {{FILE}}.html
@compile_web FILE:
    emcc {{FILE}}.c -Os -Wall -Wextra -DPLATFORM_WEB -sASSERTIONS -sALLOW_MEMORY_GROWTH -I/usr/local/include -L/usr/local/lib ./libraylib.a  -s USE_GLFW=3 -o {{FILE}}.html 
@debug FILE:
    cc {{FILE}}.c -Og -g3 -Wall -Wextra -l:libraylib.a -lm -pthread -o {{FILE}}
    just _debug {{FILE}}
