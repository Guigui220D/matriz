# Matriz

Matrix template class for linear algebra in zig.

## What is done

This is a work in progress, only simple operations are ready as of today.

- [x] Comparison
- [x] Adding/substracting
- [x] Multiplying/scaling
- [x] Transpose
- [x] Power (only for integer powers)
- [ ] Exponential of a matrix
- [ ] Rank, determinant
- [ ] Submatrices
- [ ] The rest...

## Example

```zig
var a = Mat2.init(.{
    .{1, -1},
    .{1, 0},
});

a = a.pow(4);

a = a.add(Mat2.identity);

const b = Matrix(2, 3).init(.{
    .{1, 10, 0},
    .{0, -1, 2},
});

const result = a.mul(b).scale(3);
```