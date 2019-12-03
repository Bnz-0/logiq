# Operators (or Quantum Gates)

An operator serves to modify the state of a quantum state.  
To create a new operator you have to specify a matrix which represent the operator:
```python
myOp = Op([
    [-1,0],
    [0, 1]
    ])
```
To apply an operator to a quantum state you have 2 way:
- Use the method `apply(Op)`:
    ```python
    q = qbit(1,0)
    q.apply(myOp)
    print(q) # -1|0>
    ```
- Use the `|` operator:
    ```python
    q = qbit(1,0)
    myOp | q
    print(q) # -1|0>
    ```

By default exists some common operator:
- `Op.I`
- `Op.X`
- `Op.Y`
- `Op.Z`
- `Op.H`
- `Op.cnot`
- `Op.swap`
- `Op.sqrtSwap`
- `Op.phaseGate(phase, deg)`
- `Op.C(U)` (the controlled-U gate)

---
To define an operator exists another way, that is **specify how it modify the state of a basis**.

For instance you can define an operator U in this way:
```
U such that:
    U|0> = |1>
    U|1> = |0>
```
You are surely glad to know that in logiq you can do this!
```python
U = Op.build({
    '|0>' : '|1>',
    '|1>' : '|0>'
})
```
To build an operator you have to specify a dictionary `{state : new state}`

`state` and `new state` can be expressed in different way:
- vector
- number of state of a certain basis
- string
```python
U = Op.build({
    0 : '|1>',
    1 : ket(1,0)
})
```
Of course you can specify a basis and work with it:
```python
U = Op.build({
    '|+>' : '-|->',
    '|->' : '-|+>'
}, basis=hadamard)
q = qbit('|+>', basis=hadamard)
U | q
print(q) # -1|->
```
This method can be really usefull to create dinamically an _oracle_ (see [this](Examples/Examples_list.md) for an example of use)