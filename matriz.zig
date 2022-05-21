const std = @import("std");

const T = f64;

pub fn Matrix(comptime lines: usize, comptime columns: usize) type {
    return struct {
        const Self = @This();
        const TransposeType = Matrix(columns, lines);

        const width = columns;
        const height = lines;

        comptime {
            if (lines == 0 or columns == 0)
                @compileError("Matrix must have at least one line and one column");
        }

        fn SquareMatricesMixin() type {
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

                    pub fn isInversible(self: Self) bool {
                        return (self.det() != 0);
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
                    
                    pub fn det(self: Self) T {
                        return switch (lines) {
                            0 => unreachable,
                            1 => self.det1(),
                            2 => self.det2(),
                            3 => self.det3(),
                            else => self.detN()
                        };
                    }

                    //TODO: these should be hidden: how?

                    pub fn det1(self: Self) T {
                        if (lines != 1)
                            @compileError("det1 only works on 1x1 matrices");
                        return self.data[0];
                    }

                    pub fn det2(self: Self) T {
                        if (lines != 2)
                            @compileError("det2 only works on 2x2 matrices");
                        return self.data[0] * self.data[3] - self.data[1] * self.data[2];
                    }

                    pub fn det3(self: Self) T {
                        if (lines != 3)
                            @compileError("det3 only works on 3x3 matrices");
                        var ret: T = 0;
                        ret += self.data[0] * self.data[4] * self.data[8];
                        ret += self.data[1] * self.data[5] * self.data[6];
                        ret += self.data[2] * self.data[3] * self.data[7];
                        ret -= self.data[2] * self.data[4] * self.data[6];
                        ret -= self.data[1] * self.data[3] * self.data[8];
                        ret -= self.data[0] * self.data[5] * self.data[7];
                        return ret;
                    }

                    pub fn detN(self: Self) T {
                        var i: usize = 0;
                        var ret: T = 0;
                        while (i < columns) : (i += 1) {
                            const submat = self.subMatrix(0, i);
                            var cofactor = submat.det() * self.get(0, i);
                            if (i % 2 == 1)
                                cofactor = -cofactor;
                            ret += cofactor;
                        }
                        return ret;
                    }
                };
            } else {
                return struct {};
            }
        }

        pub usingnamespace SquareMatricesMixin();

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

        pub fn sub(a: Self, b: Self) Self {
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

        // TODO: remove
        pub fn equals(a: Self, b: Self) bool {
            for (a.data) |va, i| {
                const vb = b.data[i];
                if (vb != va)
                    return false;
            }
            return true;
        }

        fn expectEqual(expected: Self, actual: Self, tolerance: T) !void {
            for (expected.data) |va, i| {
                const vb = actual.data[i];
                if (!std.math.approxEqRel(T, va, vb, tolerance)){
                    std.debug.print("Matrices differ. expected: {d}, found {d}\n", .{ expected, actual });
                    return error.TestExpectedEqual;
                }
            }
        }

        pub fn rank(self: Self) usize {
            if (isSquare and self.det() != 0)
                return lines;
            // TODO
        }

        const SubmatType = Matrix(lines - 1, columns - 1);
        pub fn subMatrix(self: Self, i: usize, j: usize) SubmatType {
            std.debug.assert(i < lines);
            std.debug.assert(j < columns);

            var ret: SubmatType = undefined;
            var k: usize = 0;
            while (k < lines) : (k += 1) {
                if (k == i)
                    continue;
                var l: usize = 0;
                while (l < columns) : (l += 1) {
                    if (l == j)
                        continue;
                    ret.set(if (k > i) k - 1 else k, if (l > j) l - 1 else l, self.get(k, l));
                }
            }
            return ret;
        }

        pub fn format(
            self: Self,
            comptime fmt: []const u8,
            options: std.fmt.FormatOptions,
            writer: anytype
        ) !void {
            _ = options;
            _ = fmt;

            try writer.writeAll("[ ");
            for (self.data) |v, i| {
                if (i != 0 and i % columns == 0)
                    try writer.writeAll("; ");
                try writer.print("{} ", .{v});
            }
            try writer.writeByte(']');
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

    try Matrix(3, 2).expectEqual(Matrix(3, 2).init(.{
        .{3, 6},
        .{3, 6},
        .{3, 7},
    }), c, 0.01);

    c = c.sub(a.transpose());

    try Matrix(3, 2).expectEqual(b, c, 0.01);
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

    try Matrix(2, 3).expectEqual(Matrix(2, 3).init(.{
        .{1,  8, 4},
        .{3, 26, 8},
    }), d, 0.01);

    const e = a.scale(2);

    try Mat2.expectEqual(Mat2.init(.{
        .{2, 4},
        .{6, 8},
    }), e, 0.01);
}

test "Matrix properties" {
    var id = Mat2.identity;
    try std.testing.expect(id.isDiagonal());
    try std.testing.expect(id.isInversible());

    id.set(0, 0, 0);
    id.set(0, 1, 42);
    try std.testing.expect(!id.isDiagonal());
    try std.testing.expect(!id.isInversible());
}

test "Matrix power" {
    const a = Mat2.init(.{
        .{1, -1},
        .{1, 0},
    });
    const b = a.pow(4);

    try Mat2.expectEqual(Mat2.init(.{
        .{-1, 1},
        .{-1, 0},
    }), b, 0.01);

    try Mat2.expectEqual(Mat2.identity, Mat2.identity.pow(10), 0.0001);
}

test "Submatrix" {
    const a = Mat3.init(.{
        .{6, 1, 1},
        .{4, -2, 5},
        .{2, 8, 7},
    });
    const b = a.subMatrix(1, 1);
    try Mat2.expectEqual(Mat2.init(.{
        .{6, 1},
        .{2, 7},
    }), b, 0.01);
}

test "Matrix determinant" {
    const a = Mat2.init(.{
        .{1, 2},
        .{3, 4},
    });
    try std.testing.expectApproxEqRel(@as(T, -2), a.det(), 0.001);

    const b = Mat3.init(.{
        .{6, 1, 1},
        .{4, -2, 5},
        .{2, 8, 7},
    });
    try std.testing.expectApproxEqRel(@as(T, -306), b.det(), 0.001);

    const c = Mat4.init(.{
        .{1, 0, 0, 0},
        .{0, 1, 3, 0},
        .{0, 0, 3, 2},
        .{0, 0, 0, 1},
    });
    _ = c;
    //try std.testing.expectApproxEqRel(@as(T, 3), c.det(), 0.001);

    const d = Mat(6).init(.{
        .{1, 0, 0, 0, 0, 0},
        .{0, 2, 3.5, 0, 2, 0},
        .{0, 2, 7, 0, 0, 0},
        .{0, 0, 3, 0.3, 0, 0},
        .{0, 0, 3, 0, 7, 0},
        .{0, 0, 0.1, -1, -1, 2},
    });
    try std.testing.expectApproxEqRel(@as(T, 36.6), d.det(), 0.001);
}