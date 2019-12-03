# Logi<span style="color:red">q</span>: "how to use"

Let's see some simple example to understand how `Qstate`, `Basis` and `Op` work togheder

> For [Qstate](Quantum_state_creations.md), [Basis](Bases.md) and [Operator](Operators.md) (Op) introduction follow the links.

---

### Start with the most used operators: the _**pauli gates**_!

Surely you know how [pauli gates works](https://en.wikipedia.org/wiki/Quantum_logic_gate#Pauli-X_gate), so go on and try to do something in logiq.
```python
# first create our qubit for testing the pauli operators
q = qbit(1,0) # |0>

# then apply the X gate
Op.X | q

# now the internal state of q has changed, but how?
# in logiq normally you can cheat and see the state without destroy it, so let's print it!
print(q)
```
The result will be `"+1|1>"` (of course), that's because the state `|1>` is represented by the vector `(0,1)`, and `X * (1,0) = (0,1)`.  
This is what was happen to the internal state of `q`.

The nice things is that **in logiq you can forget it and focus on an higher level of abstraction considering the "state in a certain basis" and not the vectors and the matrices behind this**.

> If you want to work with the matrices and vectors, you can get the internal state of a Qstate as a vector:
>```python
>q = qbit(ket(0,1)) # this is a vector so no problem
>Op.X | q
>pritn( q.state ) # in q.state there is the vector that represents the internal state of q, so the result will be "|(1+0j); 0j>"
>```
>For more on this and the interaction between logiq and numpy, please see [this](numpy_integration.md).

Continue to play with logiq and apply the other 2 gates:
```python
Op.Y | q
Op.Z | q

# now our qubit has became...
print(q) # +1j|0>

# but I'm tired because I have write the "same thing" 3 times...
# luckily in logiq you can do much more than a single operation per time!

AllInOneOperator = Op.X * Op.Y * Op.Z
# this new operator if applied to a qubit have the same effect if I applied X, Y and Z to a qubit

# so I can rewrite all things that I'vd done before in just 3 line:
q = qbit(1,0)
AllInOneOperator | q
print(q) # +1j|0>
```

---
### Let's measure the qubits

Measure a qubit is the only way to "read" the state (destroying it etc etc..), in logiq a `Basis` is needed to do this.

> Keep in mind that you can specify the basis to use every time or set it as a "default basis" for a quantum state

Measure something means 2 things:
- the internal state of the measured qubit change
- the return value is the eigenvalue associated to the autostate measured

```python
q = qbit('|0>') # the defaul basis of q will be stdbasis
result = q.measure()

# Now in result there is the eigenvalue and the internal state of q is collapsed
print(result, q) # (1+0j) +1|0>

# All of us (because we are smart people) know that this is the only possible result. That's because we had measured a qubit in state |0> with stdbasis:
print(stdbasis)
# |0>: |(1+0j); 0j>
# |1>: |0j; (1+0j)>

# and this means that if I measure q with stdbasis the probability of collapsing are:
q.printProbs()
# |0>: 100%
# |1>: 0%

# change a little bit the things:
Op.H | q

# now q is in +0.70711|0> +0.70711|1> state, and if we measuring again with stdbasis ...
q.printProbs()
# |0>: 50%
# |1>: 50%

# ... the result will be incert
```

> The value of each probabilities are in `q.prob(i [, basis])` where `i` is the i-th autostate of the basis

If you want to use a different basis it can be done or changing the "default" basis with `.stdBasis(new default basis)` method or specifying it every time:
```python
q = qbit('|0>')
q.printProbs(hadamard)
# |+>: 50%
# |->: 50%

result = q.measure(hadamard)
```

The result of a measure is fundamental but using the stdbasis you will fall into a big problem: the eigenvalue of sdbasis (the possible result after a measurement) are 1 and... 1!  
This means that you never know in which state was collapsed your qubit.

It's time to make our basis to avoid this problem:
```python
B = Basis([
    ket(1,0),
    ket(0,-1)
], 'ab')

print(B)
# |a>: |(1+0j); 0j>
# |b>: |0j; (-1+0j)>
print(B.ew) # [(1+0j), (-1+0j)]

#now we can use this basis to measure without cheating
q = qbit('|a>', basis=B)
Op.H | q
r = q.measure()

if r == B.ew[0]: print('q was collapsed into |a>')
else : print('q was collapsed into |b>')
```

---
### Adding another qubit to make things more intresting

Re-start creating 2 different qubit.  
Now, you have to know that **in logiq all qubit have a life of its own**, so the concept of _registers_ and _position of the qubits_ dows not exists here!

```python
#creating q0 and q1 with same state (|0>)
q0, q1 = qbit(1,0), qbit(1,0)

# now we can apply some gates to q0 and q1 separately ...
Op.X | q0

#... but what if I want to apply a "bigger" operator?
# To apply, for example, the operator CNOT to q0 and q1 we have to create a new quantum state:
q = q0 @ q1

# q is the "tensor product" between q0 and q1.
print(q) # "+1|10>" = |1> ⊗ |0>

#This new quantum state contains 2 qubit but we can use it as a "normal" quantum state
Op.cnot | q

# doing this we are applying the CNOT gate to q0 and q1 (using q0 as checker and q1 as target)
# now we can check the state both watching the q or q0 and q1:
print(q0) # ±1|1>
print(q1) # ±1|1>
print(q) # +1|11>
```
If you notice when q0 and q1 are printed a `±` was appeared, it means that q0 and q1 are in entanglement and what you see is an aproximated state.

> **NB:** to create 2 qubit with same state you cannot do
>```python
>q0 = q1 = qbit(1,0)
>```
>because doing this q0 **is** q1 and you have created only 1 qubit.


Remember that kind of state are `Qstate`, so all method that you can use for a simple qubit work with a much complex state.

Thanks to this, all method that you can use for a simple qubit work with a much complex state.

For example the "composition" between quantum states using the tensor product is flexile because `q0` and `q` are both `Qstate`:
```python
q = qbit(1,0) @ qbit(0,1)
print(q) # +1|01>

q = qbit(1,0) @ q
print(q) # +1|001>

q @= qbit(1,1, normalize=True)
print(q) # +0.70711|0010> +0.70711|0011>
```

> For obvious reasons "compose" two time the same quantum state (`q @ q`) cause an error

---
## The abstraction of logiq

For the last part of this fiendly introduction, I want you to understand what is the level of abstraction that I've put in logiq.

Probably the most imporant layer of abstracion was the **elimination of the concept of registers**.

In fact this help a lot to manage problem that appear in classic way, for instance in a quantum computer you are not able to apply a q-gate to 2 arbitrary qubit:
```python
# I'll take this opportunity to show you some shortcut
q = qbit('|1>') @ 2 # |11>

Op.H | q[0] # apply H only to qubit in position 0
Op.Z ^ q # apply Z to all qubit in q

print(q) # -0.70711|01> -0.70711|11>

# now the problem: imagine that we have 3 qubit...
q @= qbit(1,0)
print(q) # -0.70711|010> -0.70711|110>

#... and we want to apply a CNOT gate to the first and the last qubit.
# In a classic way can be hard or impossible do something like that, but not here:
Op.cnot | (q[0] @ q[2])

# We have create another quantum state q[0] @ q[2] and we applied to this the CNOT gate without any effort

print(q) # ±0.5|010> ±0.5|011> ±0.5|110> ±0.5|111>
```

Another abstraction is that doesn't matter what is the space of the quantum states, they work well even with different spaces:
```python
X3 = Op.build({
    0: 1,
    1: 2,
    2: 0
})

q = qbit(1,0,0) @ 3

X3 | q[1]
X3*X3 | q[2]

print(q) # +1|012>
```

For other (more complicated) example follow [this link](Examples/Examples_list.md)!