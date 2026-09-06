const std = @import("std");
const fasta = @import("root.zig");



pub fn main(init: std.process.Init) !void {
    var queue: std.Io.Queue(fasta.Fasta) = .init(&.{});
    const stdout = std.Io.File.stdout();
    
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 3) {
        std.debug.print("Usage: {s} <dna|protein> <filename>\n", .{args[0]});
        return;
    }

    var bmtype: fasta.biomolecule = undefined;
    if (std.meta.stringToEnum(fasta.biomolecule, args[1])) |value| {
        bmtype = value;
    } else {
        std.debug.print("No biomolecule type '{s}'\n", .{args[1]});
        return;
    }

    const filepath = args[2];
    const file = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer file.close(init.io);

    var producer_task = try init.io.concurrent(
        fasta.parse3,
        .{ init.io, init.gpa, &queue, file }
    );
    defer producer_task.cancel(init.io) catch {};
    var counter: u16 = 0;

    const t_start = std.Io.Timestamp.now(init.io, .awake);
    
    while (true) {
        var myFasta = queue.getOne(init.io) catch |err| switch (err) {
            error.Closed => break,
            error.Canceled => return,
        };
        defer myFasta.deinit(init.gpa);
        counter += 1;

        try stdout.writeStreamingAll(init.io, "---------\n");
        const hprint = try std.fmt.allocPrint(
            init.gpa,
            "{s:>9} {s}\n",
            .{ "header:", myFasta.header }
        );
        defer init.gpa.free(hprint);
        try stdout.writeStreamingAll(init.io, hprint);

        const sprint = try std.fmt.allocPrint(
            init.gpa,
            "sequence: {s}\n",
            .{myFasta.sequence}
        );
        defer init.gpa.free(sprint);
        try stdout.writeStreamingAll(init.io, sprint);
        
        switch (bmtype) {
            .protein => {
                const p = try fasta.Protein.init(&myFasta);

                const massprint = try std.fmt.allocPrint(
                    init.gpa,
                    "mass: {d:>0.2}\n",
                    .{p.mass}
                );
                defer init.gpa.free(massprint);
                try stdout.writeStreamingAll(init.io, massprint);
            },
            .dna => {
                const d = try fasta.DNA.init(&myFasta, init.gpa);
                defer d.deinit(init.gpa);

                const revprint = try std.fmt.allocPrint(
                    init.gpa,
                    "revcomp: {s}\n",
                    .{d.complement}
                );
                defer init.gpa.free(revprint);
                try stdout.writeStreamingAll(init.io, revprint);
            },
        }
    }

    // const elapsed: u64 =
    //     @intCast(t_start.durationTo(std.Io.Timestamp.now(init.io, .awake)).milliseconds);
    const elapsed = t_start.durationTo(std.Io.Timestamp.now(init.io, .awake)).toMilliseconds();
    const elapsedPrint = try std.fmt.allocPrint(
        init.gpa,
        "elapsed time: {d} mS\n",
        .{elapsed}
    );
    defer init.gpa.free(elapsedPrint);
    try stdout.writeStreamingAll(init.io, elapsedPrint);
    
    const finalTally = try std.fmt.allocPrint(
        init.gpa,
        "Created {d} Fasta objects\n",
        .{counter}
    );
    defer init.gpa.free(finalTally);
    try stdout.writeStreamingAll(init.io, finalTally);
}
