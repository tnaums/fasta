const std = @import("std");

const DNA = struct {
    fasta: *Fasta,
    complement: []const u8,

    pub fn init(fasta: *Fasta) !DNA {
        return .{
            .fasta = fasta,
            .complement = try DNA.reverseComplement(fasta.allocator, fasta.sequence),
        };
    }

    pub fn deinit(self: DNA) void {
        const allocator = self.fasta.allocator;
        allocator.free(self.complement);
    }

    fn reverseComplement(allocator: std.mem.Allocator, forward: []const u8) ![]const u8 {
    var revcomp = try allocator.alloc(u8, forward.len);
    var i = forward.len;
    while (i > 0) : (i -= 1) {
        var newchar: u8 = undefined;
        switch (forward[i-1]) {
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

const Protein = struct {
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

const Fasta = struct {
    header: []u8,
    sequence: []u8,
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, header: []const u8, sequence: []const u8) !Fasta {
        const h = try allocator.dupe(u8, header);
        errdefer allocator.free(h);

        const s = try allocator.dupe(u8, sequence);
        errdefer allocator.free(s);

        return .{
            .header = h,
            .sequence = s,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: Fasta) void {
        const allocator = self.allocator;
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


fn readFasta(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(Fasta), file: std.Io.File) !void {
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
                        const f: Fasta = try .init(allocator, header.items, sequence.items);
                        try queue.putOne(io, f);
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
    const f = try Fasta.init(allocator, header.items, sequence.items);
    try queue.putOne(io, f);
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

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename>\n", .{args[0]});
        return;
    }
    
    var queue: std.Io.Queue(Fasta) = .init(&.{});
    const filepath = args[1];
    const file = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer file.close(init.io);

    
    var producer_task = try init.io.concurrent(readFasta, .{ init.io, init.gpa, &queue, file });
    defer producer_task.cancel(init.io) catch {};
    var counter: u16 = 0;

    while (true) {
        var myFasta = queue.getOne(init.io) catch |err| switch (err) {
            error.Closed => break,
            error.Canceled => return,
        };
        defer myFasta.deinit();
        counter += 1;
//        const p = try Protein.init(&myFasta);
        const d = try DNA.init(&myFasta);
        defer d.deinit();
        
        try stdout.writeStreamingAll(init.io, "---------\n");
        const hprint = try std.fmt.allocPrint(init.gpa, "{s:>9} {s}\n", .{"header:", d.fasta.header});
        defer init.gpa.free(hprint);
        try stdout.writeStreamingAll(init.io, hprint);

        const sprint = try std.fmt.allocPrint(init.gpa, "sequence: {s}\n", .{d.fasta.sequence});
        defer init.gpa.free(sprint);
        try stdout.writeStreamingAll(init.io, sprint);

        const revprint = try std.fmt.allocPrint(init.gpa, "revcomp: {s}\n", .{d.complement});
        defer init.gpa.free(revprint);
        try stdout.writeStreamingAll(init.io, revprint);


    }

    const finalTally = try std.fmt.allocPrint(init.gpa, "Created {d} Fasta objects\n", .{counter});
    defer init.gpa.free(finalTally);
    try stdout.writeStreamingAll(init.io, finalTally);
}
