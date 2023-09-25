```@meta
CurrentModule = NamedRowArrays
```

# NamedRowArrays

[NamedRowArrays](https://github.com/lucasvalenzuela/NamedRowArrays.jl) provides an array type that knows about its
row names.
This allows for indexing by name in a straightforward manner in addition to the regular indexing possibilities.

These arrays are broadcastable, but note that the broadcasting behavior is undefined when broadcasting over
two equally sized NamedRowArrays with different row names. In such cases, one or both of the arrays should be
collected: `collect(A)`

```@autodocs
Modules = [NamedRowArrays]
```
