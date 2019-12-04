# Logi<span style="color:red">q</span>: quantum state creation

- [Creation of a qubit](#creation-of-a-qubit)
- [Creation of _qutrit_ and more](#creation-of-qutrit-and-more)
- [State with custom basis](#state-with-custom-basis)
- [Create a quantum state using **Qstate**](#create-a-quantum-state-using-Qstate)
- [Parsing states from string](#parsing-states-from-string)
- [Composed quantum states](#composed-quantum-states)

---
## Creation of a **qubit**
You can create a qubit with any (allowed) state you want just specify the internal state of this new qubit.  
There are mainly 2 ways to create a qubit:
- Form the vectorial representation of internal state:
    ```python
    q = qbit(ket(0,1)) 	# |1>
    q = qbit(1,0) 		# |0>
    q = qbit(0,-1j)		# -1j|1>
    ```
    > `j` is the imaginary number in python

    In the first example we have passed the ket vector to qbit function, but is equivalent to pass directly the value of the internal state.

- From a string representation:
    ```python
    q = qbit("|0>") 		# |0>
    q = qbit("-1|0> +0|1>")	# -1|0>
    ```
    This method is really useful if you want to teach something and you do not want to go into details, but it has a problem:  
    What if you want to create the state `1/√2(|0>+|1>)`?  
    You can easily build it defining the vector:
    ```python
    rd = 1/sqrt(2)
    q = qbit(rd,rd) # +0.70711|0> +0.70711|1>
    ```
    But if you want to write it as a string you have to do something more complex.

    To avoid this exists a way to say: "whatever value I insert, you normalize them and create a correct qubit"
    ```python
    q = qbit("1|0> +1|1>") # IllegalOperationError: A quantum state must have norm 1, this one has norm (2+0j)

    q = qbit("1|0> +1|1>", normalize=True) # +0.70711|0> +0.70711|1>
    ```

---
## Creation of **qutrit** and more
In logiq you are not limited to create only qubit (i.e. a **bidimensional quantum state**), in fact if for some reason, which could be physical research or implementation of a powerful communication protocol, you can create a quantum state with more than a 2-dimensional state.

The creation of this state is really simple because is identical to the [qubit creation](#creation-of-a-qubit) but just with more state:
```python
# qutrits
q = qbit(0,0,1) # |2>
q = qbit(1,2,3, normalize=True) # +0.26726|0> +0.53452|1> +0.80178|2>

# A qu... I don't know how naming this
q = qbit("0.7|6> -0.7|9>", normalize=True) # +0.70711|6> -0.70711|9>
```
Simple, not?
> **NB**: the third quantum state created is a 10-dimensional state because the 10th state (`|9>`) is the higher I have mentioned, but for a better explanation read [this](#Parsing-states-from-string)

---
## State with custom basis
In logiq every quantum state has a **basis**, it serves to measure and represent the state and by default this basis is `{|0>,|1>}`, which is called `stdbasis`.

> Or better, by default it will create a basis with the same dimension of the state.  
> So doing `qbit(0,1,0)` you create a 3-dimensional quantum state with the basis `{|0>, |1>, |2>}` associated

Of course you can create a quantum state using the basis you want but what does it mean to "_associate a basis to a state_"?


- **Measurement**

    When you measure a state without specifying the basis the "associated" basis will be used to measure
    ```python
    q = qbit(1,0, basis=hadamard)
    q.measure() # it measure using hadamard {|+>,|->}
    q.measure(stdbasis) # this time it uses the stdbasis basis to measure
    ```
- **Representation**

    When you print a quantum state it uses the associated basis to represent it, instead to forcing another representation you need to use the `printAs(basis)` method
    ```python
    q = qbit(0,1, basis=hadamard)
    print(q) # +0.70711|+> -0.70711|->
    q.printAs(stdbasis) # +1|1>
    ```
To change the associated basis, you have to use the `setBasis(basis)` method:
```python
q = qbit(1,0)
print(q) # +1|0>
q.setBasis(hadamard)
print(q) # +0.70711|+> +0.70711|->
```

---
## Create a quantum state using **Qstate**
The function `qbit` is just a shortcut to create a `Qstate` object.  
Both result is equivalent, only change the syntax to use.

- creation from vector
    ```python
    q = Qstate(ket(1,0)) #|0>, as qbit() functions
    ```
    > Notice that `Qstate(1,0)` **does not work**!
- creation from string
    ```python
    q = Qstate.parse("|1>") #as qbit() functions
    ```
> The `normalize` and the `basis` argument work in both methods, but they can be passed just using the _right position_ and not necessarily specifying it:
> ```python
> q = Qstate(ket(2,1), stdbasis, True) #+0.89443|0> +0.44721|1>
> q = Qstate.parse("|+>", hadamard, False) # |+> (normalize=False is unnecessary)
> ```


---
## Parsing states from string
When you parse a quantum state from string you have to keep in mind how it works:

### <u>Representation</u>
The representation is the same used by physics and computer scientist when talking of quantum computing: the [Bra-ket notation](https://en.wikipedia.org/wiki/Bra%E2%80%93ket_notation).
> "`a|x> + b|y> + ...`" where `a` and `b` are complex number and `x` and `y` are states of a certain basis.

### <u>Default assumptions</u>
By default, when you write a quantum state as a string, the parser assume that you are writing a state represented in the _stdbasis_: `{|0>, |1>}`  
So the state "`-1|0>`" is parsed in the state `-1|0> +0|1>`, which is probably what you want to mean

If instead you write "`0.7|0> + 0.7|3>`" the situation changes: in fact, the parser changes its assumption because `|3> ∉ {|0>, |1>}` and try to understand the written state using this different basis: `{|0>, |1>, |2>, |3>}`  
That because all symbols ("0" and "3") are in the STD_SYMBOLS.

> In the _STD_SYMBOLS_ there are **all numbers** and **all capital letters**: 0,1,...,8,9,A,B,...,Y,Z

It works in this way because if you write "`x|1> + y|4>`" you probably want to mean a 5-dimensional state where the basis is `{|0>, |1>, |2>, |3>, |4>}` and **not** a 2-dimensional state in basis `{|1>, |4>}`

> A tricky way to write n-dimensional state without specifying the basis is to write "`...` **`+0|n>`**":  
>for instance the 10-dim state `|3>` can be create with the string "`|3> +0|9>`"

### <u>Custom states</u>
You can write any kind of state without specify the basis, the behaviour is to create a new basis with exactly the symbols you use in the state.

For example, if you want to "work" with electron spin, which could have 2 possible configurations: **up** and **down**, you can easily write the state and the output will be what you expect 😉:
```python
e = qbit('1|up> +0|down>')
print(e) # +1|up>
print(e.basis)
#|up>: |(1+0j); 0j>
#|down>: |0j; (1+0j)>
```
> You must specify **all configurations** (or a basis)

**You can do this only if you are not interested in what the basis is**!  
That because the new basis generated are a simple _CanonBasis_ and you haven't any control on the eigenstate of the basis.

> Last trick: because the capital letters are in STD_SYMBOLS, if you want to generate a general state without specify basis you can simply do this:
> ```python
> q = qbit('1|a> +0|b>')
> ```
> This generates a similar state using `{|a>, |b>}` as a basis.  
> If you use the capital letter you obtain a 12-dimensional state
> ```python
> >>> len(qbit('1|a> +0|b>'))
> 2
> >>> len(qbit('1|A> +0|B>'))
> 12
> ```


---
## Composed quantum states
You can create quantum states which are the "composition" of other quantum states.

The "default" way to create this kind of states is to use the _tensor product_ (`@`) but you can also create this using the functions `qbit()` and `Qstate()`:
```python
q = qbit(
    qbit(1,0),
    Qstate(ket(0,1)),
    qbit(1,1,normalize=True)
)
print(q) # +0.70711|010> +0.70711|011>

q = Qstate( #notice the list!
    [
        qbit(1,0),
        Qstate(ket(0,1)),
        qbit(1,1,normalize=True)
    ]
)
print(q) # +0.70711|010> +0.70711|011>
```
Of course you can also compose both "_single state_" and "_multiple state_" in the same way:
```python
q = Qstate(
    [
        qbit(qbit(1,0), Qstate(ket(0,1))), #|01>
        qbit(1,1,normalize=True) #+0.70711|0> +0.70711|1>
    ]
)
print(q) #+0.70711|010> +0.70711|011>
```
But this way is not really readable, so I suggest you to use the tensor product (except for particular cases, of course)
