const Lexer = @import("lexer.zig").Lexer;

const ast = @import("ast.zig");
const std = @import("std");

const assert = std.debug.assert;

pub const Parser = struct {
    lexer: Lexer,
    nodes: std.ArrayList(ast.Node),
    allocator: std.mem.Allocator,

    //set on init
    current_token: Lexer.Token = undefined,
    peek_token: Lexer.Token = undefined,

    //begin public interface
    //create a new parser
    pub fn init(input: []const u8, allocator: std.mem.Allocator) !Parser {
        var ret = Parser{
            .lexer = try Lexer.init(input),
            .allocator = allocator,
            .nodes = std.ArrayList(ast.Node).init(allocator),
        };

        ret.nextToken();
        ret.nextToken();
        return ret;
    }

    //free parser memory
    pub fn deinit(self: *Parser) void {
        for (self.nodes.items) |node| {
            switch (node) {
                .block_expression => |b| b.expressions.deinit(),
                else => {},
            }
        }
        self.nodes.deinit();
    }

    //wipe the parser to an initial state
    pub fn reset(self: *Parser, new_input: []const u8) !void {
        self.nodes.clearRetainingCapacity();
        self.lexer = try Lexer.init(new_input);

        self.nextToken();
        self.nextToken();
    }

    //kick off parser
    pub fn parse(self: *Parser) void {
        while (self.current_token.type != .EOF) : (self.nextToken()) {
            const statement = switch (self.current_token.type) {
                .NEG => self.parseNegation(),
                else => self.parseExpression(),
            };
            _ = statement;
        }
    }

    pub fn print(self: *const Parser) void {
        for (self.nodes.items, 0..) |node, idx| {
            std.debug.print("[id={d} ", .{idx});
            self.printNode(node, 0);
            std.debug.print("]\n", .{});
        }
    }

    //end public interface
    fn printNode(self: *const Parser, node: ast.Node, indent: u32) void {
        printIndent(indent);
        switch (node) {
            .atom => |x| self.lexer.print(x.symbol, std.debug),
            .infix_expression => |x| self.printInfix(x, indent + 1),
            .block_expression => |x| {
                std.debug.print("(", .{});
                for (x.expressions.items) |n| {
                    std.debug.print("\n", .{});
                    self.printNode(self.getNode(n).*, indent + 1);
                }
                std.debug.print(")", .{});
            },
            .negation => |x| {
                std.debug.print("(", .{});
                self.lexer.print(x.symbol, std.debug);
                std.debug.print("\n", .{});
                self.printNode(self.getNode(x.right).*, indent + 1);
                std.debug.print(")", .{});
            },
            else => {},
        }
    }

    fn printIndent(indent: u32) void {
        for (0..indent) |_| std.debug.print("\t", .{});
    }

    fn printInfix(self: *const Parser, infix: ast.InfixExpression, indent: u32) void {
        const left = switch (infix) {
            inline else => |x| self.getNode(x.left),
        };
        const right = switch (infix) {
            inline else => |x| self.getNode(x.right),
        };
        const symbol = switch (infix) {
            inline else => |x| x.symbol,
        };

        std.debug.print("(", .{});
        self.lexer.print(symbol, std.debug);

        std.debug.print("\n", .{});
        self.printNode(left.*, indent);

        std.debug.print("\n", .{});
        self.printNode(right.*, indent);

        std.debug.print("\n", .{});
        printIndent(indent);
        std.debug.print(")\n", .{});
    }

    //move the parsing cursor forward by a token, setting the current and look ahead token
    fn nextToken(self: *Parser) void {
        self.current_token = self.peek_token;
        self.peek_token = self.lexer.nextToken();
    }

    //copy a given node off the stack and into heap memory via an array list
    fn addNewNode(self: *Parser, T: type, new_node: T) usize {
        const node = switch (T) {
            ast.Atom => ast.Node{ .atom = new_node },
            ast.InfixExpression => ast.Node{ .infix_expression = new_node },
            ast.Negation => ast.Node{ .negation = new_node },
            ast.BlockExpression => ast.Node{ .block_expression = new_node },
            ast.Eof => ast.Node{ .eof = new_node },
            else => {
                std.debug.print("Unhandled type: {s}\n", .{@typeName(T)});
                @panic("Unhandled type!");
            },
        };

        self.nodes.append(node) catch @panic("OOM");
        return self.nodes.items.len - 1;
    }

    //return a node from the array list based on index
    fn getNode(self: *const Parser, index: usize) *ast.Node {
        assert(index < self.nodes.items.len);

        return &self.nodes.items[index];
    }

    //given a token type return a potential prefix parsing function to handle it
    fn getPrefixFn(token: Lexer.Token.Type, peek: Lexer.Token.Type) ?*const fn (*Parser) usize {
        return switch (token) {
            //if the token after a symbol is a paren, assume its a predicate like 'P(x)' instead of just an atom 'x'
            .SYMBOL => if (peek != .LPAREN) parseAtom else null,
            .LPAREN => parseBlockExpression,
            .NEG => parseNegation,
            .EOF => parseEof,
            else => null,
        };
    }

    //given a token type return a potential infix parsing function to handle it
    fn getInfixFn(token: Lexer.Token.Type) ?*const fn (*Parser, usize) usize {
        return switch (token) {
            .IFF, .IF, .AND, .OR => parseInfixExpression,
            else => null,
        };
    }

    //entry point for parsing
    fn parseExpression(self: *Parser) ?usize {
        const prefix_fn_maybe = getPrefixFn(self.current_token.type, self.peek_token.type);
        if (prefix_fn_maybe == null) {
            std.debug.print("No prefix fn found for {s}\n", .{@tagName(self.current_token.type)});
            return null;
        }

        const prefix_fn = prefix_fn_maybe.?;
        var left_expr = prefix_fn(self);

        while (self.peek_token.type != .EOF) {
            const infix_fn_maybe = getInfixFn(self.peek_token.type);
            if (infix_fn_maybe == null) {
                return left_expr;
            }

            self.nextToken();
            const infix_fn = infix_fn_maybe.?;
            left_expr = infix_fn(self, left_expr);
        }
        return left_expr;
    }

    //expressions that has an expression to the left and right of it
    fn parseInfixExpression(self: *Parser, left: usize) usize {
        const current_token = self.current_token;
        self.nextToken();

        //set right to 0, so the infix expression's index is before whatever is right of it
        const infix_expr = switch (current_token.type) {
            .IFF => ast.InfixExpression{ .biconditional = .{ .left = left, .symbol = current_token, .right = 0 } },
            .IF => ast.InfixExpression{ .conditional = .{ .left = left, .symbol = current_token, .right = 0 } },
            .AND => ast.InfixExpression{ .conjunction = .{ .left = left, .symbol = current_token, .right = 0 } },
            .OR => ast.InfixExpression{ .disjunction = .{ .left = left, .symbol = current_token, .right = 0 } },
            else => {
                std.debug.print("Unhandled type {s}\n", .{@tagName(current_token.type)});
                @panic("");
            },
        };

        const index = self.addNewNode(ast.InfixExpression, infix_expr);
        const expr = &self.getNode(index).infix_expression;

        switch (expr.*) {
            //todo: handle null expression case
            inline .conditional, .biconditional, .conjunction, .disjunction => |*x| x.right = self.parseExpression().?,
        }

        return index;
    }

    //handle '~'
    fn parseNegation(self: *Parser) usize {
        const index = self.addNewNode(ast.Negation, ast.Negation{ .symbol = self.current_token, .right = 0 });
        const node = self.getNode(index);

        self.nextToken();
        node.negation.right = self.parseExpression().?;

        return index;
    }

    //handle symbols
    fn parseAtom(self: *Parser) usize {
        const atom = ast.Atom{ .symbol = self.current_token };
        return self.addNewNode(ast.Atom, atom);
    }

    //handle '<->'
    fn parseBiConditional(self: *Parser, left: usize) usize {
        const current_token = self.current_token;
        self.nextToken();

        const bi_cond = ast.Biconditional{ .left = left, .symbol = current_token, .right = 0 };
        const index = self.addNewNode(ast.Biconditional, bi_cond);
        const infix_expr = self.getNode(index);
        infix_expr.biconditional.right = self.parseExpression().?;

        return index;
    }

    //handle '->'
    fn parseConditional(self: *Parser, left: usize) usize {
        const current_token = self.current_token;
        self.nextToken();

        //set right to 0, so the infix expression's index is before whatever is right of it
        const cond = ast.Conditional{ .left = left, .symbol = current_token, .right = 0 };
        const index = self.addNewNode(ast.Conditional, cond);
        const infix_expr = self.getNode(index);
        infix_expr.conditional.right = self.parseExpression().?;

        return index;
    }

    //handle parens '( ... )'
    fn parseBlockExpression(self: *Parser) usize {
        //consume '('
        self.nextToken();
        const block = ast.BlockExpression{ .expressions = std.ArrayList(usize).init(self.allocator) };
        const index = self.addNewNode(ast.BlockExpression, block);
        var node = self.getNode(index);

        while (self.current_token.type != .RPAREN and self.current_token.type != .EOF) : (self.nextToken()) {
            const expr_maybe = self.parseExpression();
            if (expr_maybe == null) continue;

            const expr = expr_maybe.?;
            node.block_expression.expressions.append(expr) catch @panic("OOM");
        }

        //todo handle err
        if (self.current_token.type != .RPAREN) @panic("Expected )");
        return index;
    }

    fn parseEof(self: *Parser) usize {
        const last = ast.Eof{ .symbol = self.current_token };
        return self.addNewNode(ast.Eof, last);
    }

    fn parsePredicate(_: *Parser) void {}
};

//tests

test "negation" {
    var parser = try Parser.init("~x", std.testing.allocator);
    defer parser.deinit();
    parser.parse();

    try std.testing.expectEqual(2, parser.nodes.items.len);
    const atom_1 = parser.nodes.items[1].atom;
    const neg = parser.nodes.items[0].negation;

    try std.testing.expectEqualStrings("x", atom_1.symbol.getLiteral(parser.lexer.input));
    try std.testing.expectEqualDeep(atom_1, parser.getNode(neg.right).atom);

    try parser.reset("~(x -> y)");
    parser.parse();

    const block = parser.getNode(neg.right).block_expression;
    try std.testing.expectEqual(1, block.expressions.items.len);

    const infix = parser.getNode(block.expressions.items[0]).infix_expression.conditional;
    try std.testing.expectEqualStrings("x", parser.getNode(infix.left).atom.symbol.getLiteral(parser.lexer.input));
    try std.testing.expectEqualStrings("y", parser.getNode(infix.right).atom.symbol.getLiteral(parser.lexer.input));
}

test "conditionals" {
    var parser = try Parser.init("x -> y", std.testing.allocator);
    defer parser.deinit();
    parser.parse();

    try std.testing.expectEqual(3, parser.nodes.items.len);
    const atom_1 = parser.nodes.items[0].atom;
    const atom_2 = parser.nodes.items[2].atom;

    try std.testing.expectEqualStrings("x", atom_1.symbol.getLiteral(parser.lexer.input));
    try std.testing.expectEqualStrings("y", atom_2.symbol.getLiteral(parser.lexer.input));

    const cond = parser.nodes.items[1].infix_expression.conditional;
    try std.testing.expectEqualDeep(atom_1, parser.getNode(cond.left).atom);
    try std.testing.expectEqualDeep(atom_2, parser.getNode(cond.right).atom);

    try parser.reset("a <-> b -> c");
    parser.parse();
    try std.testing.expectEqual(5, parser.nodes.items.len);

    const atom_3 = parser.nodes.items[0].atom;
    const atom_4 = parser.nodes.items[2].atom;
    try std.testing.expectEqualStrings("a", atom_3.symbol.getLiteral(parser.lexer.input));
    try std.testing.expectEqualStrings("b", atom_4.symbol.getLiteral(parser.lexer.input));

    const bi_cond = parser.nodes.items[1].infix_expression.biconditional;
    try std.testing.expectEqualDeep(atom_3, parser.getNode(bi_cond.left).atom);

    const cond_2 = parser.getNode(bi_cond.right).infix_expression.conditional;
    try std.testing.expectEqualDeep(atom_4, parser.getNode(cond_2.left).atom);

    try parser.reset("(x ^ y) -> z");
    parser.parse();
    parser.print();
    try std.testing.expectEqual(6, parser.nodes.items.len);
}

test "and / or" {
    var parser = try Parser.init("x ^ y v x", std.testing.allocator);
    defer parser.deinit();
    parser.parse();

    const conj = parser.nodes.items[1].infix_expression.conjunction;
    try std.testing.expectEqualStrings("x", parser.getNode(conj.left).atom.symbol.getLiteral(parser.lexer.input));
}

test "printing" {
    var parser = try Parser.init("(x ^ y) v x", std.testing.allocator);
    defer parser.deinit();
    parser.parse();
    parser.print();
}
