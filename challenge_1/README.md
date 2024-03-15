# How to run
First at all: 
```bash
cd \path\to\this\folder; cd challenge_1
```

All the object files (.o) will be store in the obj/ directory. 
To compile all use: 
```bash
make
``` 

To configure it use the file config.hpp where you can setup all the parameters, the function and the gradient to use. There is the possibility to change the return type of the function and the gradient, do it carefully.

If you change some parameters in _config.hpp_, I recommend to use this: 
```bash
make clean; make
```

And after that, to run the executable: 
```bash
./main
```

# Parameters 

Strategy parameters:

- __mode:__ Define which strategy use to update $x_k$, already implemented are: _Heavy-Ball_, _Nesterov_ and _Default_ which is [Default mode formula](./img/challenge_1/Default_mode.png)
- __grad_mode:__ Define if use the gradient define or approximate it
- __strategy:__ Define which strategy use to update $\alpha_k$, already implemented are: _exponential decay_, _Inverse decay_, _Approximate Line Search (Armijo rule)_; 

Stopping criteria:

- __max_iter:__ n° of maximum iteration;
- __step_length:__ minimum distance between x_k and $x_{k+1}$;
- __residual:__ minimum norm for the gradient of $f(x)$;

Parameters needed by some mode or strategy: 

- __mu:__ parameter for the exponential decay strategy; 
- __sigma:__ parameter for the Armijo rule in the Approximate Line Search strategy; 
- __nu:__ parameter for the Heavy-ball mode; 

Initial parameters: 

- __alpha_0:__ starting $\alpha_k$;
- __x0:__ starting $x_k$;

Parameter for the approximate gradient:

- __h:__ Step-lentgh for the approximate gradient

# Code guidelines

In the __gradient_descent__ function in the __main.cpp__ the structure of the code is generally: 

```
Variables' declaration and definition;

If necessary, define x_1;

while(not reach a satisfactory result and iteration are above the max)

    update the stopping variable;

    if(satisfactory result)
        exit
    else
        update x

```