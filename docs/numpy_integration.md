# Integration with `numpy`

Every logiq's object have the `.npm()` method that returns a `numpy.matrix` object.

```python
>>> M = Op.Y.npm()
>>> type(M)
numpy.matrix

>>> stdbasis.npm() == matrix.Id(2).npm()
matrix([[ True,  True],
        [ True,  True]])
```


> This matrix is the **real** data of the object so <u>**do not modify it!**</u>  
>You can easily copy it if you need.

The **only** object with a different behaviour is `Qstate`.

To get the state as a `numpy.matrix` you have to get the internal state as a `vector` and then call the `.npm()`:
```python
>>> q = Qstate.random(2)
>>> intState = q.state #this contains the internal state as a vector
>>> intState.npm()
matrix([[-0.41113669+0.59933783j],
        [ 0.62821559-0.27767959j]])
```

The mainly motivations of this are to avoid involuntary modifications of the internal state of a `Qstate`:

```python
>>> q = Qstate.random(2)
>>> intState = q.state # this is a copy of the real internal state
>>> Op.Y | q
>>> intState.npm() == q.state.npm()
matrix([[False],
        [False]])
```

---
## Creation using numpy's objects

Except for `Qstate` every other classes generate their data using `numpy.matrix(input)`, so every type of input accepted from numpy.matrix could be a good input.

The classes `vector`, `matrix`, `Op` and `Basis` initialize their content in a way similar to this:
```python
class obj:
    def __init__(self, input):
        self.content = numpy.matrix(input)
```
So every input allowed from `numpy.matrix` are possible good input for the class and this cause a great integration with numpy's classes!