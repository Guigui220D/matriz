const std = @import("std");

const T = f64;

pub fn Matrix(comptime lines: usize, comptime columns: usize) type {
    return struct {
        const Self = @This();
        const TransposeType = Matrix(columns, lines);

        pub const height = lines;
        pub const width = columns;

        fn SquareMatricesClosure() type {
            if (isSquare) {
                return struct {
                    pub const identity = blk: {
                        var id = std.mem.zeroes(Self);
                        var i: usize = 0;
                        while (i < lines) : (i += 1) {
                            id.set(i, i, 1);
                        }
                        break :blk id;
                    };

                    pub fn isDiagonal(self: Self) bool {
                        var i: usize = 0;
                        while (i < lines) : (i += 1) {
                            var j: usize = 0;
                            while (j < columns) : (j += 1) {
                                if (j == i)
                                    continue;
                                if (self.get(i, j) != 0)
                                    return false;
                            }
                        }
                        return true;
                    }

                    // TODO: pow for non integer values
                    pub fn pow(self: Self, n: usize) Self {
                        if (n == 0) return identity;
                        if (n == 1) return self;
                        var ret = self;
                        if (self.isDiagonal()) {
                            var j: usize = 0;
                            const p = @intToFloat(T, n);
                            while (j < columns) : (j += 1) {
                                ret.set(j, j, std.math.pow(T, self.get(j, j), p));
                            }
                        } else {
                            // Terribly expensive for big values of n
                            // TODO: other method
                            var i: usize = 1;
                            while (i < n) : (i += 1) {
                                ret = ret.mul(self);
                            }
                        }
                        return ret;
                    }

                    // TODO: exp
                };
            } else {
                return struct {};
            }
        }

        pub usingnamespace SquareMatricesClosure();

        const isSquare = (lines == columns);
        
        pub fn transpose(self: Self) TransposeType {
            var ret: TransposeType = undefined;
            var i: usize = 0;
            while (i < lines) : (i += 1) {
                var j: usize = 0;
                while (j < columns) : (j += 1) {
                    ret.set(j, i, self.get(i, j));
                }
            }
            return ret;
        }

        /// Create a new matrix
        /// Example for Mat2:
        /// Mat2.init(.{ .{1, 2}, .{3, 4} })
        pub fn init(contents: anytype) Self {
            var ret: Self = undefined;

            comptime var max_line: usize = 0;
            inline for (contents) |line, i| {
                if (i >= lines)
                    @compileError("Too many lines in matrix");
                comptime var max_column: usize = 0;
                inline for (line) |elem, j| {
                    if (j >= columns)
                        @compileError("Too many columns in matrix");
                    ret.set(i, j, elem);
                    max_column = j;
                }
                if (max_column < columns - 1)
                    @compileError("Too few columns in matrix");
                max_line = i;
            }
            if (max_line < lines - 1)
                @compileError("Too few lines in matrix");

            return ret;
        }

        pub fn set(self: *Self, i: usize, j: usize, value: T) void {
            self.data[i + j * lines] = value;
        }

        pub fn get(self: Self, i: usize, j: usize) T {
            return self.data[i + j * lines];
        }

        pub fn add(a: Self, b: Self) Self {
            var ret: Self = undefined;
            for (a.data) |va, i| {
                const vb = b.data[i];
                ret.data[i] = va + vb;
            }
            return ret;
        }

        pub fn sub (a: Self, b: Self) Self {
            var ret: Self = undefined;
            for (a.data) |va, i| {
                const vb = b.data[i];
                ret.data[i] = va - vb;
            }
            return ret;
        }

        fn MulResType(comptime other: type) type {
            if (columns != other.height)
                @compileError("Incompatible matrix dimensions");
            return Matrix(lines, other.width);
        }

        pub fn mul(a: Self, b: anytype) MulResType(@TypeOf(b)) {
            var ret = std.mem.zeroes(MulResType(@TypeOf(b)));
            var i: usize = 0;
            while (i < lines) : (i += 1) {
                var j: usize = 0;
                while (j < @TypeOf(b).width) : (j += 1) {
                    var k: usize = 0;
                    var sum: T = 0;
                    while (k < columns) : (k += 1) {
                        sum += a.get(i, k) * b.get(k, j);
                    }
                    ret.set(i, j, sum);
                }
            }
            return ret;
        }

        pub fn scale(a: Self, b: T) Self {
            var ret: Self = undefined;
            for (a.data) |va, i| {
                ret.data[i] = va * b;
            }
            return ret;
        }

        pub fn equals(a: Self, b: Self) bool {
            for (a.data) |va, i| {
                const vb = b.data[i];
                if (vb != va)
                    return false;
            }
            return true;
        }

        data: [lines * columns]T,
    };
}

pub const Mat2 = Matrix(2, 2);
pub const mat2 = Mat2.init;
pub const Mat3 = Matrix(3, 3);
pub const mat3 = Mat3.init;
pub const Mat4 = Matrix(4, 4);
pub const mat4 = Mat4.init;

pub fn Mat(comptime size: usize) type {
    return Matrix(size, size);
}

test "Creating a matrix, adding, transposing" {
    const a = Matrix(2, 3).init(.{
        .{1, 2, 3},
        .{4, 5, 6},
    });

    const b = Matrix(3, 2).init(.{
        .{2, 2},
        .{1, 1},
        .{0, 1},
    });

    var c = b.add(a.transpose());

    try std.testing.expect(c.equals(Matrix(3, 2).init(.{
        .{3, 6},
        .{3, 6},
        .{3, 7},
    })));

    c = c.sub(a.transpose());

    try std.testing.expect(c.equals(b));
}

test "Multiplication" {
    const id = Mat2.identity;
    const a = Mat2.init(.{
        .{1, 2},
        .{3, 4},
    });

    const b = a.mul(id);

    try std.testing.expect(b.equals(a));

    const c = Matrix(2, 3).init(.{
        .{1, 10, 0},
        .{0, -1, 2},
    });

    const d = a.mul(c);

    try std.testing.expect(d.equals(Matrix(2, 3).init(.{
        .{1,  8, 4},
        .{3, 26, 8},
    })));

    const e = a.scale(2);
    try std.testing.expect(e.equals(Matrix(2, 2).init(.{
        .{2, 4},
        .{6, 8},
    })));
}

test "Matrix properties" {
    var id = Matrix(5, 5).identity;
    try std.testing.expect(id.isDiagonal());

    id.set(1, 2, 42.0);
    try std.testing.expect(!id.isDiagonal());
}

test "Matrix power" {
    const a = Mat2.init(.{
        .{1, -1},
        .{1, 0},
    });
    const b = a.pow(4);
    try std.testing.expect(b.equals(Mat2.init(.{
        .{-1, 1},
        .{-1, 0},
    })));

    try std.testing.expect(Mat2.identity.equals(Mat2.identity.pow(10)));
}