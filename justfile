_run FILE:
    ./{{FILE}}
_debug FILE:
    valgrind ./{{FILE}}
@run-fast FILE:
    cc {{FILE}}.c -O3 -Wall -Wextra -pedantic -ffast-math -funsafe-math-optimizations -l:libraylib.a -lm -pthread -o {{FILE}}
    just _run {{FILE}}
@run FILE:
    cc {{FILE}}.c -O3 -Wall -Wextra -pedantic -l:libraylib.a -lm -pthread -o {{FILE}}
    just _run {{FILE}}
@debug FILE:
    cc {{FILE}}.c -Og -g3 -Wall -Wextra -pedantic -l:libraylib.a -lm -pthread -o {{FILE}}
    just _debug {{FILE}}
