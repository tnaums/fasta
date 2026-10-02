//! By convention, root.zig is the root source file when making a package.
const std = @import("std");
const Io = std.Io;

/// This is a documentation comment to explain the `printAnotherMessage` function below.
///
/// Accepting an `Io.Writer` instance is a handy way to write reusable code.
pub fn printAnotherMessage(writer: *Io.Writer) Io.Writer.Error!void {
    try writer.print("Run `zig build test` to run the tests.\n", .{});
}

pub const FastaError = error{
    ParsingProblem,
    InitializingFailure,
    TranslationError,
};

const geneticCode = std.StaticStringMap(u8).initComptime([_]struct { []const u8, u8 }{
    .{ "TTT", 'F' }, .{ "TTC", 'F' }, .{ "TTG", 'L' }, .{ "TTA", 'L' },
    .{ "CTT", 'L' }, .{ "CTC", 'L' }, .{ "CTA", 'L' }, .{ "CTG", 'L' },
    .{ "ATT", 'I' }, .{ "ATC", 'I' }, .{ "ATA", 'I' }, .{ "ATG", 'M' },
    .{ "GTT", 'V' }, .{ "GTC", 'V' }, .{ "GTA", 'V' }, .{ "GTG", 'V' },
    .{ "TCT", 'S' }, .{ "TCC", 'S' }, .{ "TCA", 'S' }, .{ "TCG", 'S' },
    .{ "CCT", 'P' }, .{ "CCC", 'P' }, .{ "CCA", 'P' }, .{ "CCG", 'P' },
    .{ "ACT", 'T' }, .{ "ACC", 'T' }, .{ "ACA", 'T' }, .{ "ACG", 'T' },
    .{ "GCT", 'A' }, .{ "GCC", 'A' }, .{ "GCA", 'A' }, .{ "GCG", 'A' },
    .{ "TAT", 'Y' }, .{ "TAC", 'Y' }, .{ "TAA", '*' }, .{ "TAG", '*' },
    .{ "CAT", 'H' }, .{ "CAC", 'H' }, .{ "CAA", 'Q' }, .{ "CAG", 'Q' },
    .{ "AAT", 'N' }, .{ "AAC", 'N' }, .{ "AAA", 'K' }, .{ "AAG", 'K' },
    .{ "GAT", 'D' }, .{ "GAC", 'D' }, .{ "GAA", 'E' }, .{ "GAG", 'E' },
    .{ "TGT", 'C' }, .{ "TGC", 'C' }, .{ "TGA", '*' }, .{ "TGG", 'W' },
    .{ "CGT", 'R' }, .{ "CGC", 'R' }, .{ "CGA", 'R' }, .{ "CGG", 'R' },
    .{ "AGT", 'S' }, .{ "AGC", 'S' }, .{ "AGA", 'R' }, .{ "AGG", 'R' },
    .{ "GGT", 'G' }, .{ "GGC", 'G' }, .{ "GGA", 'G' }, .{ "GGG", 'G' },
});

pub const biomolecule = enum { protein, dna };
pub const readingframe = enum { first, second, third };

pub const DNA = struct {
    header: []const u8,
    sequence: []const u8,
    complement: []const u8,
    translation: ?[6][]u8 = null,

    
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
        if (self.translation) |value| {
            for (0..6) |t| {allocator.free(value[t]);}
        }
        allocator.free(self.complement);
        allocator.free(self.sequence);
        allocator.free(self.header);
    }

    pub fn ligate(first: *DNA, second: *DNA, allocator: std.mem.Allocator) !DNA {
        const combined_header = try std.fmt.allocPrint(allocator, "{s}_+_{s}", .{ first.header, second.header });
        defer allocator.free(combined_header);
        const combined_sequence = try std.fmt.allocPrint(allocator, "{s}{s}", .{ first.sequence, second.sequence });
        defer allocator.free(combined_sequence);

        return DNA.init(allocator, combined_header, combined_sequence);
    }

    pub fn mapDNA(self: *DNA, allocator: std.mem.Allocator, io: Io, file: Io.File, frame: readingframe) !void {
        const frameInt = @intFromEnum(frame);
        const ruler = "----:----|----:----|----:----|----:----|----:----|----:----|";
        if (self.translation) |value| {
            _ = value;
        } else {
            self.translation = try self.translate(allocator);
        }
        var selectedFrame: []u8 = undefined;
        if (self.translation) |value| {
            selectedFrame = value[frameInt]; // selects the []u8
        } else {
            unreachable;
        }
        var aminoacids = std.ArrayList(u8).empty;
        defer aminoacids.deinit(allocator);
        switch(frame) {
            .first => {
                for (selectedFrame) |aa| {
                    try aminoacids.append(allocator, aa);
                    try aminoacids.append(allocator, ' ');
                    try aminoacids.append(allocator, ' ');                    
                }
            },
            .second => {
                for (selectedFrame) |aa| {
                    try aminoacids.append(allocator, ' ');
                    try aminoacids.append(allocator, aa);
                    try aminoacids.append(allocator, ' ');
                }
            },
            .third => {
                for (selectedFrame) |aa| {
                    try aminoacids.append(allocator, ' ');
                    try aminoacids.append(allocator, ' ');
                    try aminoacids.append(allocator, aa);
                }  
            },
        }
        var index: usize = 0;
        var buffer1: [70]u8 = undefined;
        var buffer2: [70]u8 = undefined;
        var buffer3: [70]u8 = undefined;
        var buffer4: [70]u8 = undefined;        
        while (index < self.sequence.len - 60) : (index += 60) {

            const line1 = try std.fmt.bufPrint(&buffer1, "{s}\n", .{self.sequence[index..index + 60]});
            const line2 = try std.fmt.bufPrint(&buffer2, "{s} {d:>5}\n", .{ ruler, index + 60 });
            const line3 = try std.fmt.bufPrint(&buffer3, "{s}\n", .{self.complement[index..index + 60]});
            const line4 = try std.fmt.bufPrint(&buffer4, "{s}\n", .{aminoacids.items[index..index + 60]});
            try file.writeStreamingAll(io, line1);
            try file.writeStreamingAll(io, line2);
            try file.writeStreamingAll(io, line3);
            try file.writeStreamingAll(io, line4);
            try file.writeStreamingAll(io, "\n");
        }
        const line1 = try std.fmt.bufPrint(&buffer1, "{s}\n", .{self.sequence[index..]});
        const line2 = try std.fmt.bufPrint(&buffer2, "{s} {d:>5}\n", .{ruler[0..line1.len - 1], index + line1.len - 1});
        const line3 = try std.fmt.bufPrint(&buffer3, "{s}\n", .{self.complement[index..]});
        const line4 = try std.fmt.bufPrint(&buffer4, "{s}\n", .{aminoacids.items[index..]});
        try file.writeStreamingAll(io, line1);
        try file.writeStreamingAll(io, line2);
        try file.writeStreamingAll(io, line3);
        try file.writeStreamingAll(io, line4);
    }

    pub fn addTranslation(self: *DNA, allocator: std.mem.Allocator) !void {
        self.translation = try self.translate(allocator);
    }

    pub fn printOrfs(self: DNA, allocator: std.mem.Allocator, cutoff: ?usize) void {
        _ = allocator;
        var filter: usize = undefined;
        if (cutoff) |value| {
            filter = value;
        } else {
            filter = 0;
        }
        
        var seq: []u8 = undefined;
        if (self.translation) |value| {
            for (0..6) |i| {
                seq = value[i];
                var orfs = std.mem.tokenizeScalar(u8, seq, '*');
                while (true) {
                    if (orfs.next()) |orf| {
                        if (orf.len > filter) {
                            std.debug.print("frame {d}: ", .{i});
                            std.debug.print("{s}\n\n", .{orf});
                        }
                    } else {
                        break;
                    }
                }
            }
        }

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

    pub fn format(self: DNA, writer: *Io.Writer) !void {
        const bp = self.sequence.len;
        try writer.print(">{s}|{d}bp\n", .{ self.header, bp });
        var lineIndex: usize = 0;
        while (lineIndex + 60 < self.sequence.len) : (lineIndex += 60) {
            try writer.print("{s}\n", .{self.sequence[lineIndex .. lineIndex + 60]});
        }
        try writer.print("{s}\n", .{self.sequence[lineIndex..]});
    }

    fn translate(self: DNA, allocator: std.mem.Allocator) ![6][]u8 {
        var accumulator: [6][]u8 = undefined;
        var j: usize = 0;

        while (j < 3) : (j += 1) {
            var aminoacids = std.ArrayList(u8).empty;
            defer aminoacids.deinit(allocator);

            var index: usize = 0 + j;
            while (index < self.sequence.len - 2) : (index += 3) {
                const codon = self.sequence[index .. index + 3];
                const aminoacid = geneticCode.get(codon) orelse 'X';
                try aminoacids.append(allocator, aminoacid);
            }
            accumulator[j] = try aminoacids.toOwnedSlice(allocator);
        }
        while (j < 6) : (j += 1) {
            var aminoacids = std.ArrayList(u8).empty;
            defer aminoacids.deinit(allocator);

            var index: usize = j - 3;
            while (index < self.complement.len - 2) : (index += 3) {
                const codon = self.complement[index .. index + 3];
                const aminoacid = geneticCode.get(codon) orelse 'X';
                try aminoacids.append(allocator, aminoacid);
            }
            accumulator[j] = try aminoacids.toOwnedSlice(allocator);
        }
        return accumulator;
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
        pub fn format(self: Protein, writer: *Io.Writer) !void {
            try writer.print(">{s}|{d:.2}kDa\n", .{ self.header, self.mass });
            var lineIndex: usize = 0;
            while (lineIndex + 60 < self.sequence.len) : (lineIndex += 60) {
                try writer.print("{s}\n", .{self.sequence[lineIndex .. lineIndex + 60]});
            }
            try writer.print("{s}\n", .{self.sequence[lineIndex..]});
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

/// Creates DNA structs from fasta formatted sequence files.
///
/// This function parses fasta files, creates DNA structs, and places them
/// into a queue.
///
/// # Parameters
/// - `io`: Io instance for async
/// - `queue`: - Io.Queue(DNA)
/// - `file`: - Io.File with fasta formatted DNA sequences
///
/// # Example
/// ```zig
/// var queue: std.Io.Queue(fasta.DNA) = .init(&.{});
/// var producer_taskNew = try init.io.concurrent(fasta.parseDNA, .{ init.io, init.gpa, &queue, file });
/// defer producer_taskNew.cancel(init.io) catch {};
/// while (true) {
///     var myDNA = queue.getOne(init.io) catch |err| switch (err) {
///     error.Closed => break,
///     error.Canceled => return,
///     };
///}
///```
pub fn parseDNA(io: Io, allocator: std.mem.Allocator, queue: *Io.Queue(DNA), file: Io.File) !void {
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

/// Creates Protein structs from fasta formatted sequence files.
///
/// This function parses fasta files, creates Protein structs, and places them
/// into a queue.
///
/// # Parameters
/// - `io`: Io instance for async
/// - `queue`: - Io.Queue(Protein)
/// - `file`: - Io.File with fasta formatted DNA sequences
///
/// # Example
/// ```zig
/// var queue: std.Io.Queue(fasta.Protein) = .init(&.{});
/// var producer_task = try init.io.concurrent(fasta.parseProtein, .{ init.io, init.gpa, &queue, file });
/// defer producer_task.cancel(init.io) catch {};
/// while (true) {
///     var myProtein = queue.getOne(init.io) catch |err| switch (err) {
///     error.Closed => break,
///     error.Canceled => return,
///     };
///     defer myProtein.deinit(init.gpa);
///}
///```
pub fn parseProtein(io: Io, allocator: std.mem.Allocator, queue: *Io.Queue(Protein), file: Io.File) !void {
    const state = enum { inHeader, inSequence };
    var myState: ?state = null;

    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);

    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    defer queue.close(io);
    var buf: [64]u8 = undefined;
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

pub fn parseProteinSIMD(io: std.Io, allocator: std.mem.Allocator, queue: *std.Io.Queue(Protein), file: Io.File) !void {
    defer queue.close(io);
    const state = enum { inHeader, inSequence };
    var myState: ?state = null;

    var vcounter: u32 = 0;

    // 64KB buffer for high throughput
    var buffer: [64 * 1024]u8 = undefined;
    // SIMD Configuration
    const VSize = std.simd.suggestVectorLength(u8).?;
    const Vector = @Vector(VSize, u8);
    const nl: Vector = @splat('\n');
    const ns: Vector = @splat('>');

    var header: std.ArrayList(u8) = .empty;
    var sequence: std.ArrayList(u8) = .empty;
    defer header.deinit(allocator);
    defer sequence.deinit(allocator);
    const Bitmask = std.meta.Int(.unsigned, VSize);

    while (true) {
        const n = file.readStreaming(io, &.{&buffer}) catch |err| {
            if (err == error.EndOfStream) break;
            return err;
        };
//        std.debug.print("n = {d}\n", .{n});

        var i: usize = 0;

        // 1. SIMD LOOP
        while (i + VSize <= n) : (i += VSize) {
            vcounter += 1;

            const v: Vector = buffer[i..][0..VSize].*;
            const newLines = v == nl;
            const newSequences = v == ns;

            const mask: Bitmask = @bitCast(newSequences);
            const mask2: Bitmask = @bitCast(newLines);

            const newsequences_index_left = @ctz(mask); // matches left-most '>'
            const newlines_index_left = @ctz(mask2); // matches left-most '\n'
            const newlines_number: u8 = @popCount(mask2); // might be 0 or 1 or 2 or 3 new lines

            // For debugging
//            std.debug.print("\nLoop number {d}; {any}\n", .{ vcounter, myState });
//            std.debug.print("{b:>32}: {d}\n", .{ mask, newsequences_index_left });
//            std.debug.print("{b:>32}: {d}\n", .{ mask2, newlines_index_left });
            
            if (myState) |s| {
                switch (s) {
                    .inHeader => {
                        if (newlines_index_left != 0) {
                            try header.appendSlice(allocator, buffer[i .. i + newlines_index_left]);
                        }
                        try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + VSize]);
                        myState = state.inSequence;
                    },
                    .inSequence => {
                        if (newsequences_index_left != VSize) { // found a new sequence
                            if (newlines_number == 3) {
                                const newlines_index_right = @clz(mask2); // matches right-most '\n'
                                try sequence.appendSlice(allocator, buffer[i .. i + newlines_index_left]);
                                try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + newsequences_index_left - 1]);
                                // put sequence in queue here
                                const p = try Protein.init(allocator, header.items, sequence.items);
                                //                            defer p.deinit(allocator);
                                try queue.putOne(io, p);

                                header.clearRetainingCapacity();
                                sequence.clearRetainingCapacity();
                                try header.appendSlice(allocator, buffer[i + newsequences_index_left + 1 .. i + VSize - 1 - newlines_index_right]);
                                try sequence.appendSlice(allocator, buffer[i + VSize - newlines_index_right .. i + VSize]);
                                continue;
                            }
                            if (newlines_number == 2) {
                                if (newlines_index_left + 1 > newsequences_index_left - 1) {
                                    const newlines_index_right = @clz(mask2); // matches right-most '\n'
                                    try sequence.appendSlice(allocator, buffer[i .. i + newlines_index_left]);

                                    const p = try Protein.init(allocator, header.items, sequence.items);
                                    try queue.putOne(io, p);

                                    header.clearRetainingCapacity();
                                    sequence.clearRetainingCapacity();

                                    try header.appendSlice(allocator, buffer[i + newsequences_index_left + 1 .. i + VSize - 1 - newlines_index_right]);
                                    try sequence.appendSlice(allocator, buffer[i + VSize - newlines_index_right .. i + VSize]);
                                    continue;
                                }
                                try sequence.appendSlice(allocator, buffer[i .. i + newlines_index_left]);
                                try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + newsequences_index_left - 1]);
                                const p = try Protein.init(allocator, header.items, sequence.items);
                                try queue.putOne(io, p);

                                header.clearRetainingCapacity();
                                sequence.clearRetainingCapacity();

                                try header.appendSlice(allocator, buffer[i + newsequences_index_left + 1 .. i + VSize]);
                                myState = state.inHeader;
                                continue;
                            } // end of nelines_number == 2

                            if (newsequences_index_left == 0) { // starts with >
                                // put protein in queue first
                                const p = try Protein.init(allocator, header.items, sequence.items);
                                try queue.putOne(io, p);

                                header.clearRetainingCapacity();
                                sequence.clearRetainingCapacity();

                                try header.appendSlice(allocator, buffer[i + 1 .. i + newlines_index_left]);
                                try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + VSize]);
                                continue;
                            }
                            // start with sequence and ends in header; sequence\n>startOfHeader
                            try sequence.appendSlice(allocator, buffer[i .. i + newsequences_index_left - 1]);
                            const p = try Protein.init(allocator, header.items, sequence.items);
                            try queue.putOne(io, p);

                            header.clearRetainingCapacity();
                            sequence.clearRetainingCapacity();

                            try header.appendSlice(allocator, buffer[i + newsequences_index_left + 1 .. i + VSize]);
                            myState = state.inHeader;
                            continue;
                        }
                        // if there aren't any new sequence characters
                        if (newlines_number == 2) {
                            // only happens if the final character is a new line
                            try sequence.appendSlice(allocator, buffer[i .. i + newlines_index_left]);
                            try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + VSize - 1]);
                            continue;
                            
                        }
                        if (newlines_index_left == VSize) {
                            // in sequence and no newlines; just append to sequence
                            try sequence.appendSlice(allocator, buffer[i .. i + VSize]);
                            continue;
                        }
                        if (newlines_index_left == 0) {
                            // first character is \n, just append the rest
                            try sequence.appendSlice(allocator, buffer[i + 1 .. i + VSize]);
                            continue;
                        }
                        // single internal newline, append either side to sequence
                        try sequence.appendSlice(allocator, buffer[i .. i + newlines_index_left]);
                        try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + VSize]);
                    },
                }
            } else {
                if (newsequences_index_left != VSize) { // There is a new sequence
                    if (newlines_index_left == VSize) { // newsequence without newline
                        try header.appendSlice(allocator, buffer[i .. i + VSize]);
                        myState = state.inHeader;
                        continue;
                    }
                    try header.appendSlice(allocator, buffer[i + newsequences_index_left + 1 .. i + newlines_index_left]);
                    try sequence.appendSlice(allocator, buffer[i + newlines_index_left + 1 .. i + VSize]);
                    myState = state.inSequence;
                }
            }
        } // end of SIMD loop
        // process the remaining bytes in the buffer
        while (i < n) : (i += 1) {
            const byte: u8 = buffer[i];
            if (myState) |s| {
                switch (s) {
                    .inHeader => {
                        if (byte == '\n') {
                            myState = state.inSequence;
                            continue;
                        }
                        try header.append(allocator, byte);
                    },
                    .inSequence => {
                        if (byte == '>') {
                            const f: Protein = try .init(allocator, header.items, sequence.items);
                            try queue.putOne(io, f);
                            sequence.clearRetainingCapacity();
                            header.clearRetainingCapacity();
                            myState = state.inHeader;
                            continue;
                        }
                        if (byte != '\n') {
                            try sequence.append(allocator, byte);
                        }
                    },
                }
            } else {
                if (byte == '>') {
                    myState = state.inHeader;
                    continue;
                }
            }
        }
    }
    const f = try Protein.init(allocator, header.items, sequence.items);
    try queue.putOne(io, f);
}



pub fn parseSingleDNA(io: Io, allocator: std.mem.Allocator, filepath: []const u8) !DNA {
    // open file
    const file = try std.Io.Dir.cwd().openFile(io, filepath, .{});
    defer file.close(io);

    const state = enum { inHeader, inSequence };
    var myState: ?state = null;
    // initialize ArrayLists
    var header = std.ArrayList(u8).empty;
    defer header.deinit(allocator);
    var sequence = std.ArrayList(u8).empty;
    defer sequence.deinit(allocator);

    var buf: [1024]u8 = undefined;
    // parse out the info
    while (true) {
        const n = file.readStreaming(io, &.{&buf}) catch |err| {
            if (err == error.EndOfStream) break;
            return err;
        };
        if (n == 0) {
            if (header.items.len == 0) break;
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
                        if (buf[i] == '>') { // if a second sequence begins, just return the first one
                            return try DNA.init(allocator, header.items, sequence.items);
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
    return try DNA.init(allocator, header.items, sequence.items);
}

test "create Protein" {
    const testing = std.testing;
    const header = "fveg_042069";
    const sequence = "ACDEFGHIRRSTYWPNMNMYILC";
    const p = try Protein.init(std.testing.allocator, header, sequence);
    defer p.deinit(std.testing.allocator);
    try testing.expectEqualStrings(p.header, "fveg_042069");
    try testing.expect(p.mass > 2.81 and p.mass < 2.83);
    try testing.expect(p.sequence.len == 23);
}

test "create DNA" {
    const testing = std.testing;
    const header = "fveg_042069";
    const sequence = "TACTACTATTGCCAGCATTGCTGCTAAAGAAGAAGGGGTATCTCTCGAGAAAAGAGAGGCTGAAGCTCACCACCATCATCATCACCACCACGAGAATTTATACTTTCAAGCTCCTGCCGA";
    const d = try DNA.init(std.testing.allocator, header, sequence);
    defer d.deinit(std.testing.allocator);
    try testing.expectEqualStrings(d.header, "fveg_042069");
    try testing.expect(d.complement.len == 120);
    try testing.expect(std.mem.startsWith(u8, d.complement, "TCGGCAGGAG"));
    try testing.expect(std.mem.endsWith(u8, d.complement, "AATAGTAGTA"));
}

test "ligate DNA" {
    const testing = std.testing;
    const header = "fveg_042069";
    const sequence = "TACTACTATTGCCAGCATTGCTGCTAAAGAAGAAGGGGTATCTCTCGAGAAAAGAGAGGCTGAAGCTCACCACCATCATCATCACCACCACGAGAATTTATACTTTCAAGCTCCTGCCGA";
    var d = try DNA.init(testing.allocator, header, sequence);
    defer d.deinit(testing.allocator);
    try testing.expectEqualStrings(d.header, "fveg_042069");
    try testing.expect(d.complement.len == 120);
    try testing.expect(std.mem.startsWith(u8, d.complement, "TCGGCAGGAG"));
    try testing.expect(std.mem.endsWith(u8, d.complement, "AATAGTAGTA"));

    const header2 = "second_part";
    const sequence2 = "TACTACTATTGCCAGCATTGCTGCTAAAGAAGAAGGGGTATCTCTCGAGAAAAGAGAGGCTGAACTCACCACCATCATCATCACCACCACGAGAATTTATACTTTCAAGCTCCTGCCGA";
    var d2 = try DNA.init(testing.allocator, header2, sequence2);
    defer d2.deinit(std.testing.allocator);
    try testing.expectEqualStrings(d2.header, "second_part");
    try testing.expect(d2.complement.len == 119);
    try testing.expect(std.mem.startsWith(u8, d.complement, "TCGGCAGGAG"));
    try testing.expect(std.mem.endsWith(u8, d.complement, "AATAGTAGTA"));

    var combined = try DNA.ligate(&d, &d2, testing.allocator);
    defer combined.deinit(testing.allocator);
    try testing.expect(combined.sequence.len == 239);
}
