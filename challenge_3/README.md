## Class Mesh
I'm using a vector since it should be modified many times so it is faster than modified a map with key rows and value vectors.
By this way, the complexity of access and modification follow O(1) instead of O(log(n)).

Why not to use a vector of vectors?
Since, using a larger vector, I search inside a single vector with the correct index and I'm avoiding the double search of 2 vectors also if each has complexity O(1).

## muParser
To use this program you need to have installed muParser. This is the [guide](https://github.com/beltoforion/muparser/blob/master/Install.txt) on how to install it.