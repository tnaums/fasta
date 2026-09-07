const std = @import("std");
const fasta = @import("root.zig");

pub fn main(init: std.process.Init) !void {
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

    var counter: u16 = 0;

    // Starting protein queue
    var queueNew: std.Io.Queue(fasta.Protein) = .init(&.{});

    const fileNew = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer fileNew.close(init.io);

    var producer_taskNew = try init.io.concurrent(fasta.parseProtein, .{ init.io, init.gpa, &queueNew, fileNew });
    defer producer_taskNew.cancel(init.io) catch {};

    const t_start = std.Io.Timestamp.now(init.io, .awake);
    while (true) {
        var myProtein = queueNew.getOne(init.io) catch |err| switch (err) {
            error.Closed => break,
            error.Canceled => return,
        };
        defer myProtein.deinit(init.gpa);
        counter += 1;
        
        try stdout.writeStreamingAll(init.io, "---------\n");
        const hdprint = try std.fmt.allocPrint(init.gpa, "{s:>9} {s}\n", .{ "header:", myProtein.header });
        defer init.gpa.free(hdprint);
        try stdout.writeStreamingAll(init.io, hdprint);

        const sqprint = try std.fmt.allocPrint(init.gpa, "sequence: {s}\n", .{myProtein.sequence});
        defer init.gpa.free(sqprint);
        try stdout.writeStreamingAll(init.io, sqprint);

        const massNewprint = try std.fmt.allocPrint(init.gpa, "mass: {d:>0.2}\n", .{myProtein.mass});
        defer init.gpa.free(massNewprint);
        try stdout.writeStreamingAll(init.io, massNewprint);
    }
    const elapsed = t_start.durationTo(std.Io.Timestamp.now(init.io, .awake)).toMilliseconds();
    const elapsedPrint = try std.fmt.allocPrint(init.gpa, "elapsed time: {d} mS\n", .{elapsed});
    defer init.gpa.free(elapsedPrint);
    try stdout.writeStreamingAll(init.io, elapsedPrint);

    const finalTally = try std.fmt.allocPrint(init.gpa, "Created {d} Fasta objects\n", .{counter});
    defer init.gpa.free(finalTally);
    try stdout.writeStreamingAll(init.io, finalTally);
}
