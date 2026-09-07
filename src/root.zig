//! By convention, root.zig is the root source file when making a package.
const std = @import("std");
const Io = std.Io;

/// This is a documentation comment to explain the `printAnotherMessage` function below.
///
/// Accepting an `Io.Writer` instance is a handy way to write reusable code.
pub fn printAnotherMessage(writer: *Io.Writer) Io.Writer.Error!void {
    try writer.print("Run `zig build test` to run the tests.\n", .{});
}

pub fn add(a: i32, b: i32) i32 {
    return a + b;
}

test "basic add functionality" {
    try std.testing.expect(add(3, 7) == 10);
}

pub const biomolecule = enum { protein, dna };
pub const Biomolecule = union(biomolecule) {
    protein: Protein2,
    dna: DNA2,

    fn init(
        allocator: std.mem.Allocator,
        header: []u8,
        sequence: []u8,
        dorp: biomolecule,
    ) !Biomolecule {
        switch (dorp) {
            .protein => {
                return Biomolecule{ .protein = Protein2.init(allocator, header, sequence) };
            },
            .dna => {

            }
        }
    }

    fn deinit(
        self: Biomolecule,
        allocator: std.mem.Allocator,
    ) void {
        
    }
};

pub const DNA = struct {
    fasta: *Fasta,
    complement: []const u8,

    pub fn init(fasta: *Fasta, allocator: std.mem.Allocator) !DNA {
        return .{
            .fasta = fasta,
            .complement = try DNA.reverseComplement(allocator, fasta.sequence),
        };
    }

    pub fn deinit(self: DNA, allocator: std.mem.Allocator) void {
        allocator.free(self.complement);
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

pub const Protein2 = struct {
    fasta: *Fasta,
    mass: f32,

    pub fn init(allocator: std.mem.Allocator, header: []const u8, sequence: []const u8) !Protein {
        const fasta = try Fasta.init(allocator, header, sequence);
        return .{
            .fasta = &fasta,
            .mass = calculateMass(fasta.sequence),
        };
    }

    pub fn deinit(self: Protein2, allocator: std.mem.Allocator) void {
        self.fasta.deinit(allocator);
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

pub const Protein = struct {
    fasta: *Fasta,
    mass: f32,

    pub fn init(fasta: *Fasta) !Protein {
        return .{
            .fasta = fasta,
            .mass = calculateMass(fasta.sequence),
        };
    }

    // pub fn deinit(self: Protein) void {
    //     const allocator = self.fasta.allocator;
    //     allocator.destroy(self);
    // }

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

pub const Fasta = struct {
    header: []u8,
    sequence: []u8,

    pub fn init(allocator: std.mem.Allocator, header: []const u8, sequence: []const u8) !Fasta {
        const h = try allocator.dupe(u8, header);
        errdefer allocator.free(h);

        const s = try allocator.dupe(u8, sequence);
        errdefer allocator.free(s);

        return .{
            .header = h,
            .sequence = s,
        };
    }

    pub fn deinit(self: Fasta, allocator: std.mem.Allocator) void {
        allocator.free(self.sequence);
        allocator.free(self.header);
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


// Uses a small buffer; is much faster than using a single byte buffer, but after about 8 bytes, it seems to
// stay the same. Not slower than using parse2 which reads lines with takeDelimiter. How does
// file.readStreaming differ from reader.interface.takeDelimiter('\n')?
pub fn parse3(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(Fasta), file: std.Io.File) !void {
    const state = enum { inHeader, inSequence };
    var myState: ?state = null;

    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);

    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    defer queue.close(io);

    while (true) {
        var buf: [16]u8 = undefined;
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
                            const f: Fasta = try .init(allocator, header.items, sequence.items);
                            try queue.putOne(io, f);
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
    const f = try Fasta.init(allocator, header.items, sequence.items);
    try queue.putOne(io, f);
}

pub fn parse2(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(Fasta), file: std.Io.File) !void {
    var header: std.ArrayList(u8) = .empty;
    defer header.deinit(allocator);
    var sequence: std.ArrayList(u8) = .empty;
    defer sequence.deinit(allocator);
    var startFlag: bool = true;
    defer queue.close(io);

    var buf: [100]u8 = undefined; // must be big enough for longest line
    var reader: std.Io.File.Reader = file.reader(io, &buf);

    // Fasta parser, putting each Fasta struct into the queue
    while (try reader.interface.takeDelimiter('\n')) |line| {
        if (line.len == 0) {
            continue;
        }
        if (line[0] == '>') {
            if (!startFlag) {
                const f: Fasta = try .init(allocator, header.items, sequence.items);
                try queue.putOne(io, f);
                sequence.clearRetainingCapacity();
                header.clearRetainingCapacity();
            }
            try header.appendSlice(allocator, line);
            startFlag = false;
        } else {
            try sequence.appendSlice(allocator, line);
        }
    }

    const f: Fasta = try .init(allocator, header.items, sequence.items);
    try queue.putOne(io, f);
}

pub fn parse(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(Biomolecule), file: std.Io.File, dorp: biomolecule) !void {
    const state = enum { inHeader, inSequence };
    var myState: ?state = null;

    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);

    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    defer queue.close(io);

    while (true) {
        var byte: [1]u8 = undefined;
        const amt = file.readStreaming(io, &.{&byte}) catch |err| {
            if (err == error.EndOfStream) break;
            return err;
        };
        if (amt == 0) {
            if (header.items.len == 0) return;
            break;
        }

        if (myState) |s| {
            switch (s) {
                .inHeader => {
                    if (byte[0] == '\n') {
                        myState = state.inSequence;
                        continue;
                    }
                    try header.append(allocator, byte[0]);
                },
                .inSequence => {
                    if (byte[0] == '>') {
                        var f: Fasta = try .init(allocator, header.items, sequence.items);
                        var forqueue: Biomolecule = undefined;
                        switch (dorp) {
                            .dna => {
                                forqueue = Biomolecule{ .dna = try DNA.init(&f, allocator) };
                            },
                            .protein => {
                                forqueue = Biomolecule{ .protein = try Protein.init(&f) };                                
                            },
                        }
                        try queue.putOne(io, forqueue);
                        sequence.clearRetainingCapacity();
                        header.clearRetainingCapacity();
                        myState = state.inHeader;
                        continue;
                    }
                    if (byte[0] != '\n') {
                        //                        const c = std.ascii.toUpper(byte[0]);
                        try sequence.append(allocator, std.ascii.toUpper(byte[0]));
                    }
                },
            }
        } else {
            if (byte[0] == '>') {
                myState = state.inHeader;
                continue;
            }
        }
    }
    var f = try Fasta.init(allocator, header.items, sequence.items);
    var forqueue: Biomolecule = undefined;
    switch (dorp) {
        .dna => {
            forqueue = Biomolecule{ .dna = try DNA.init(&f, allocator) };
        },
        .protein => {
            forqueue = Biomolecule{ .protein = try Protein.init(&f) };                                
        },
    }
    try queue.putOne(io, forqueue);    
}

pub fn fastaConsumer(
    io: std.Io,
    queue: *std.Io.Queue(Fasta),
) !Fasta {
    const value = queue.getOne(io) catch |err| {
        return err;
    };
    return value;
}
