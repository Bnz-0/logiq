# Schrodinger's cat explanation

> Have you ever met this [cat](https://en.wikipedia.org/wiki/Schr%C3%B6dinger%27s_cat)?

```python
import math
from logiq import *

#first we need a basis to using as "point of view" ...
schrodinger = CanonBasis(2, ['Alive','Dead'])

print(schrodinger)
#|Alive>: |(1+0j); 0j>
#|Dead>: |0j; (1+0j)>

#... and a box that changes the state of our cat:
rd = 1/math.sqrt(2)
box = Op.build({
    #if the state is |Alive> this box "decrease" the "aliveness" and "increase" the "deadness"
    '|Alive>': '{}|Alive> {}|Dead>'.format(rd,rd),
    #for some complex reason this box have to change also the |Death> state increasing the aliveness
    '|Dead>': '{}|Alive> {}|Dead>'.format(rd,-rd)
    },
    basis=schrodinger #(just to remember that we are speaking in schrodinger language!)
)


#after that we have to find an alive cat...
cat = qbit('|Alive>', basis=schrodinger)

#... and put it inside the box
box | cat

#now the poor cat is inside the box and we didn't know if it's alive or dead, to see it we must open the box and observe...
cat.measure()

print(cat) #will our cat survive?
```
---

In **logiq** we can cheat and see the state of a qubit.  
If we do it after putting the cat in the box we will see
```python
. . .
box | cat
print(cat)
# "+0.70711|Alive> +0.70711|Dead>"
```
This means that **the cat is both alive and dead** and only when we open the box, i.e. measure, the cat will remain alive or die.
