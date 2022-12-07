### Базовые действия

#### Создание целевой функции для минимизации
```

from qubovert import boolean_var

N = 10

# create the variables
x = {i: boolean_var('x(%d)' % i) for i in range(N)}

# minimize \sum_{i=0}^{N-2} (1-2x_{i}) x_{i+1}
model = 0
for i in range(N-1):
    model += (1 - 2 * x[i]) * x[i+1]

# subject to the constraint that x_1 equals the XOR of x_3 and x_5
# enforce with a penalty factor of 3
model.add_constraint_eq_XOR(x[1], x[3], x[5], lam=3)

# subject to the constraints that the sum of all variables is less than 4
# enforce with a penalty factor of 5
model.add_constraint_lt_zero(sum(x.values()) - 4, lam=5)

```

