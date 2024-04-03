const std = @import("std");

const Parser = @import("parser.zig").Parser;

pub fn main() !void {
    const in = std.io.getStdIn();
    const out = std.io.getStdOut();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    while (true) {
        try out.writer().writeAll(">> ");
        var msg_buf: [4096]u8 = undefined;
        const msg = try in.reader().readUntilDelimiterOrEof(&msg_buf, '\n');

        if (msg) |m| {
            var parser = try Parser.init(m, allocator);
            defer parser.deinit();

            parser.parse();
        } else {
            break;
        }
    }
}
