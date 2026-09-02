const std = @import("std");
const fasta = @import("root.zig");

pub fn main(init: std.process.Init) !void {
    const stdout = std.Io.File.stdout();
    const args = try init.minimal.args.toSlice(init.arena.allocator());
    if (args.len != 2) {
        std.debug.print("Usage: {s} <filename>\n", .{args[0]});
        return;
    }
    
    var queue: std.Io.Queue(fasta.Fasta) = .init(&.{});
    const filepath = args[1];
    const file = try std.Io.Dir.cwd().openFile(init.io, filepath, .{});
    defer file.close(init.io);

    
    var producer_task = try init.io.concurrent(fasta.readFasta, .{ init.io, init.gpa, &queue, file });
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
        const d = try fasta.DNA.init(&myFasta);
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
