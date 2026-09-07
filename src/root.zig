//! By convention, root.zig is the root source file when making a package.
const std = @import("std");
const Io = std.Io;

/// This is a documentation comment to explain the `printAnotherMessage` function below.
///
/// Accepting an `Io.Writer` instance is a handy way to write reusable code.
pub fn printAnotherMessage(writer: *Io.Writer) Io.Writer.Error!void {
    try writer.print("Run `zig build test` to run the tests.\n", .{});
}

pub const biomolecule = enum { protein, dna };

pub const DNA = struct {
    header: []const u8,
    sequence: []const u8,
    complement: []const u8,

    pub fn init(allocator: std.mem.Allocator, header: []const u8, sequence: []const u8) !DNA {
        const h = try allocator.dupe(u8, header);
        errdefer allocator.free(h);

        const s = try allocator.dupe(u8, sequence);
        errdefer allocator.free(s);

        return .{
            .header = h,
            .sequence = s,
            .complement = try DNA.reverseComplement(allocator, sequence),
        };
    }

    pub fn deinit(self: DNA, allocator: std.mem.Allocator) void {
        allocator.free(self.complement);
        allocator.free(self.sequence);
        allocator.free(self.header);
    }

    fn reverseComplement(allocator: std.mem.Allocator, forward: []const u8) ![]const u8 {
        var revcomp = try allocator.alloc(u8, forward.len);
        var i = forward.len;
        while (i > 0) : (i -= 1) {
            var newchar: u8 = undefined;
            switch (forward[i - 1]) {
                'A' => newchar = 'T',
                'C' => newchar = 'G',
                'G' => newchar = 'C',
                'T' => newchar = 'A',
                else => newchar = '?',
            }
            revcomp[forward.len - i] = newchar;
        }
        return revcomp;
    }
};

pub const Protein = struct {
    header: []const u8,
    sequence: []const u8,
    mass: f32,

    pub fn init(allocator: std.mem.Allocator, header: []const u8, sequence: []const u8) !Protein {
        const h = try allocator.dupe(u8, header);
        errdefer allocator.free(h);

        const s = try allocator.dupe(u8, sequence);
        errdefer allocator.free(s);

        return .{
            .header = h,
            .sequence = s,
            .mass = calculateMass(sequence),
        };
    }

    pub fn deinit(self: Protein, allocator: std.mem.Allocator) void {
        allocator.free(self.sequence);
        allocator.free(self.header);
    }

    fn calculateMass(sequence: []const u8) f32 {
        var mass: f32 = 18.0;
        for (sequence) |aa| {
            const k = std.meta.stringToEnum(AminoAcid, &[_]u8{aa});
            if (k) |key| {
                mass += massMap.get(key);
            } else if (aa == '*') {
                return mass / 1000; // stop codon, we are done
            } else {
                return 0.0; // something went wrong
            }
        }

        return mass / 1000;
    }    
};

const AminoAcid = enum {
    A,
    C,
    D,
    E,
    F,
    G,
    H,
    I,
    K,
    L,
    M,
    N,
    P,
    Q,
    R,
    S,
    T,
    V,
    W,
    Y,
};

const massMap: std.EnumArray(AminoAcid, f32) = .init(.{
    .A = 71.07855,
    .C = 103.14464,
    .D = 115.08826,
    .E = 129.11504,
    .F = 147.17571,
    .G = 57.05177,
    .H = 137.14062,
    .I = 113.15890,
    .K = 128.17358,
    .L = 113.15890,
    .M = 131.19820,
    .N = 114.10354,
    .P = 97.11623,
    .Q = 128.13032,
    .R = 156.18707,
    .S = 87.07796,
    .T = 101.10474,
    .V = 99.13211,
    .W = 186.21220,
    .Y = 163.17512,
});


pub fn parseDNA(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(DNA), file: std.Io.File) !void {
    const state = enum { inHeader, inSequence };
    var myState: ?state = null;

    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);

    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    defer queue.close(io);
    var buf: [1024]u8 = undefined;
    while (true) {
        const n = file.readStreaming(io, &.{&buf}) catch |err| {
            if (err == error.EndOfStream) break;
            return err;
        };
        if (n == 0) {
            if (header.items.len == 0) return;
            break;
        }
        var i: u16 = 0;
        while (i < n) : (i += 1) {
            if (myState) |s| {
                switch (s) {
                    .inHeader => {
                        if (buf[i] == '\n') {
                            myState = state.inSequence;
                            continue;
                        }
                        try header.append(allocator, buf[i]);
                    },
                    .inSequence => {
                        if (buf[i] == '>') {
                            const d: DNA = try .init(allocator, header.items, sequence.items);
                            try queue.putOne(io, d);
                            sequence.clearRetainingCapacity();
                            header.clearRetainingCapacity();
                            myState = state.inHeader;
                            continue;
                        }
                        if (buf[i] != '\n') {
                            try sequence.append(allocator, std.ascii.toUpper(buf[i]));
                        }
                    },
                }
            } else {
                if (buf[i] == '>') {
                    myState = state.inHeader;
                    continue;
                }
            }
        }
    }
    const d = try DNA.init(allocator, header.items, sequence.items);
    try queue.putOne(io, d);
}

pub fn parseProtein(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(Protein), file: std.Io.File) !void {
    const state = enum { inHeader, inSequence };
    var myState: ?state = null;

    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);

    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    defer queue.close(io);
    var buf: [1024]u8 = undefined;
    while (true) {
        const n = file.readStreaming(io, &.{&buf}) catch |err| {
            if (err == error.EndOfStream) break;
            return err;
        };
        if (n == 0) {
            if (header.items.len == 0) return;
            break;
        }
        var i: u16 = 0;
        while (i < n) : (i += 1) {
            if (myState) |s| {
                switch (s) {
                    .inHeader => {
                        if (buf[i] == '\n') {
                            myState = state.inSequence;
                            continue;
                        }
                        try header.append(allocator, buf[i]);
                    },
                    .inSequence => {
                        if (buf[i] == '>') {
                            const p: Protein = try .init(allocator, header.items, sequence.items);
                            try queue.putOne(io, p);
                            sequence.clearRetainingCapacity();
                            header.clearRetainingCapacity();
                            myState = state.inHeader;
                            continue;
                        }
                        if (buf[i] != '\n') {
                            try sequence.append(allocator, std.ascii.toUpper(buf[i]));
                        }
                    },
                }
            } else {
                if (buf[i] == '>') {
                    myState = state.inHeader;
                    continue;
                }
            }
        }
    }
    const p = try Protein.init(allocator, header.items, sequence.items);
    try queue.putOne(io, p);
}

// pub fn fastaConsumer(
//     io: std.Io,
//     queue: *std.Io.Queue(Fasta),
// ) !Fasta {
//     const value = queue.getOne(io) catch |err| {
//         return err;
//     };
//     return value;
// }
