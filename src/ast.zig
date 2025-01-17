const Token = @import("lexer.zig").Lexer.Token;

const std = @import("std");

//possible types a node can be
pub const Node = union(enum) {
    atom: Atom,
    predicate: Predicate,
    negation: Negation,
    infix_expression: InfixExpression,
    block_expression: BlockExpression,
    quantifier: Quantifier,
    eof: Eof,
};

//wrapper for infix expressions (->, <->, etc)
pub const InfixExpression = union(enum) {
    conditional: Conditional,
    biconditional: Biconditional,
    conjunction: Conjunction,
    disjunction: Disjunction,
};

//a block of expressions within parens, e.g: (expr)
pub const BlockExpression = struct {
    expressions: std.ArrayList(usize),
};

//an atom is a literal, something that fits in between connectives (e.g x,y,P)
pub const Atom = struct {
    symbol: Token,
};

//marker for eof
pub const Eof = struct {
    symbol: Token, //for now use for index at end of file
};

//a predicate follows a quantifier (e.g P(x), P(x,y))
pub const Predicate = struct {
    //the variable literal right before the parens
    symbol: Token,
    //the list of args
    variables: std.ArrayList(Atom),
};

// (expr ^ expr)
pub const Conjunction = struct {
    left: usize, //index to node
    right: usize, //index to node

    symbol: Token, //not sure if this is needed
};

// (expr v expr)
pub const Disjunction = struct {
    left: usize, //index to node
    right: usize, //index to node

    symbol: Token, //not sure if this is needed
};

// ~expr
pub const Negation = struct {
    right: usize, //expression right of neg symbol

    symbol: Token, //might not be needed?
};

// expr <-> expr
pub const Biconditional = struct {
    left: usize, //expression left of <->
    right: usize, //expression right of <->

    symbol: Token, //might not be needed?
};

// expr -> expr
pub const Conditional = struct {
    left: usize, //expression left of ->
    right: usize, //expression right of ->

    symbol: Token, //might not be needed?
};

// In ast.zig, add the Quantifier node type:
pub const Quantifier = struct {
    symbol: Token, // The quantifier token (∀ or ∃)
    variable: Token, // The bound variable
    expression: usize, // The expression that follows
};
