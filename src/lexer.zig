const std = @import("std");
const assert = std.debug.assert;

pub const Lexer = struct {
    //possible errors
    const LexerErr = error{
        EOF,
        InputToLarge,
    };

    pub const Token = struct {
        //what type of token this is
        type: Type,
        //starting index on the original input source of when this token starts
        start_pos: u32,
        //when to stop reading the original input for this token end
        end_pos: u32,

        //possible token types
        pub const Type = enum {
            //metadata
            EOF,

            //a variable
            SYMBOL,

            //connectives/primatives
            IFF,
            IF,
            AND,
            OR,
            NEG,

            //quantifiers
            FORALL,
            EXISTS,

            //parenthesis/commas
            LPAREN,
            RPAREN,
            COMMA,
        };

        //lookup table to convert strings to a token type
        const token_map = std.ComptimeStringMap(Type, .{
            .{ "<->", .IFF },
            .{ "->", .IF },
            .{ "\\exists", .EXISTS },
            .{ "\\forall", .FORALL },
        });

        pub fn getLiteral(self: Token, input: []const u8) []const u8 {
            assert(self.end_pos >= self.start_pos);
            assert(self.end_pos <= input.len);

            return input[self.start_pos..self.end_pos];
        }
    };

    //src input to lex
    input: []const u8,
    //position of cursor in input
    index: u32 = 0,
    //char at current cursor position
    current_char: u8,
    //char at position 1 ahead of current_char
    peek_char: u8 = undefined,
    //is at eof
    is_eof: bool = false,

    //begin public interface

    //create a new lexer, input size must fit in a u32
    pub fn init(input: []const u8) error{InputToLarge}!Lexer {
        assert(input.len > 0);
        //max input length needs to fit in u32
        const input_test: u32 = @truncate(input.len);
        if (input_test != input.len) return LexerErr.InputToLarge;

        return Lexer{ .input = input, .current_char = input[0] };
    }

    //move cursor forwad, consuming next token as current token and setting look ahead token as well
    pub fn nextToken(self: *Lexer) Token {
        self.skipwhitespace() catch {
            self.is_eof = true;
        };

        if (self.is_eof) return Token{ .start_pos = @truncate(self.input.len), .end_pos = @truncate(self.input.len), .type = .EOF };
        const t_type = switch (self.current_char) {
            '(' => self.buildToken(.LPAREN),
            ')' => self.buildToken(.RPAREN),
            '~' => self.buildToken(.NEG),
            ',' => self.buildToken(.COMMA),
            '^' => self.buildToken(.AND),
            'v' => self.buildToken(.OR),

            else => {
                const start_index = self.index;
                const stop_index: u32 = self.readUntilDelimiter() catch blk: {
                    self.is_eof = true;
                    break :blk @truncate(self.input.len);
                };

                assert(start_index < stop_index);
                assert(stop_index <= self.input.len);

                const type_maybe = Token.token_map.get(self.input[start_index..stop_index]);
                const t_type = if (type_maybe) |t| t else .SYMBOL;

                //return here since we already read up to the next whitespace
                return Token{ .type = t_type, .start_pos = start_index, .end_pos = stop_index };
            },
        };

        self.readChar() catch |err| switch (err) {
            error.EOF => self.is_eof = true,
        };

        return t_type;
    }

    pub fn print(self: Lexer, token: Token, writer: anytype) void {
        const input = token.getLiteral(self.input);
        writer.print("{s}: [{d}, {d}] '{s}'", .{ @tagName(token.type), token.start_pos, token.end_pos, input });
    }

    //shorthand for consuming a token thats a single char
    fn buildToken(self: *const Lexer, t_type: Token.Type) Token {
        return Token{
            .type = t_type,
            .start_pos = self.index,
            .end_pos = self.index + 1,
        };
    }

    // end public interface

    //read a single char and set peek token, if at eof file returns err
    fn readChar(self: *Lexer) error{EOF}!void {
        self.index += 1;
        if (self.index >= self.input.len) return LexerErr.EOF;

        self.current_char = self.input[self.index];
        self.peek_char = if (self.index + 1 < self.input.len) self.input[self.index + 1] else 0;
    }

    //check if given char is considered whitespace
    fn isWhiteSpace(char: u8) bool {
        return (char == ' ' or char == '\t' or char == '\n' or char == '\t' or char == '\r' or char == 0);
    }

    //keep reading until the current char is not a whitespace char
    fn skipwhitespace(self: *Lexer) error{EOF}!void {
        while (isWhiteSpace(self.current_char)) : (try self.readChar()) {}
    }

    //could the given char be a valid token char
    fn isValidChar(char: u8) bool {
        //return std.ascii.isAlphanumeric(char) or char == '>' or char == '<' or char == '~' or char == '-' or char == '\\';
        //todo figure out better way of checking if a char is allowed or not
        return char != ' ' and char != '\n' and char != '\r' and char != '\t' and char != ',' and char != ')' and char != '(';
    }

    //keep reading until we hit the end of a token literal text
    fn readUntilDelimiter(self: *Lexer) error{EOF}!u32 {
        while (isValidChar(self.current_char)) : (try self.readChar()) {}
        return self.index;
    }
};

test "lexing" {
    const input = "x y -> <-> \\forall x (x,y) P() \\exists";
    var lexer = try Lexer.init(input);
    const expected = [_]Lexer.Token.Type{ .SYMBOL, .SYMBOL, .IF, .IFF, .FORALL, .SYMBOL, .LPAREN, .SYMBOL, .COMMA, .SYMBOL, .RPAREN, .SYMBOL, .LPAREN, .RPAREN, .EXISTS, .EOF };
    const literals = [_][]const u8{ "x", "y", "->", "<->", "\\forall", "x", "(", "x", ",", "y", ")", "P", "(", ")", "\\exists", "" };

    for (expected, literals) |expected_token, literal| {
        const token = lexer.nextToken();
        //lexer.print(token, std.debug);

        try std.testing.expectEqual(expected_token, token.type);
        try std.testing.expectEqualStrings(literal, input[token.start_pos..token.end_pos]);
    }
}

test "ands" {
    const input = "x ^ y v x";
    var lexer = try Lexer.init(input);
    const expected = [_]Lexer.Token.Type{ .SYMBOL, .AND, .SYMBOL, .OR, .SYMBOL, .EOF };
    const literals = [_][]const u8{ "x", "^", "y", "v", "x", "" };

    for (expected, literals) |expected_token, literal| {
        const token = lexer.nextToken();

        try std.testing.expectEqual(expected_token, token.type);
        try std.testing.expectEqualStrings(literal, token.getLiteral(input));
    }
}

test "edge cases" {
    const input = "x + y";
    var lexer = try Lexer.init(input);
    const expected = [_]Lexer.Token.Type{ .SYMBOL, .SYMBOL, .SYMBOL, .EOF };
    const literals = [_][]const u8{ "x", "+", "y", "" };

    for (expected, literals) |expected_token, literal| {
        const token = lexer.nextToken();

        try std.testing.expectEqual(expected_token, token.type);
        try std.testing.expectEqualStrings(literal, token.getLiteral(input));
    }
}
