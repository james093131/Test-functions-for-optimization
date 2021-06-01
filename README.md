Some of Optimization Test Problems

This file provide several Optimization Test Problems written in C.

1. Ackley (-40~40)
2. Rastrigin (-5.12~5.12)
3. Rosenbrock (-5~10)
4. Sphere (0~PI)
5. Michalewicz (-5.12~5.12)
6. BentCigar (-100~100)
7. Schaffer F7 (-100~100)
8. Zakharov (-100~100)
9. Griewank (-600~600)
10. Schwefel (-500~500)
11. Sum Squares (-10~10)
12. Powell (-4~5)


In order to call the 'TestingFunction.h' function in C++ , first you need to import it 
```
#include "TestingFunction.h"
```

The following lines is just creating the arguments necessary to demonstrate how to call the function.

```
# x : Solutions
x = [0,0,0,0,0,0,0,0,0,0]
# f : Fitness Value
f = [0.0]
# DIM : Dimensions (In general we will set Dimension to 2,10,30,50,100)
DIM = 10
# N : Number of Solutions
N = 1
# func_num : Function number
func_num = 1
```

Call Function 
```
Testing_Function(x, f, DIM, N, func_num)
```

