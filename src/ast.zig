const Token = @import("lexer.zig").Lexer.Token;

const std = @import("std");

//possible types a node can be
pub const Node = union(enum) {
    atom: Atom,
    predicate: Predicate,
    negation: Negation,
    infix_expression: InfixExpression,
    block_expression: BlockExpression,
};

//wrapper for infix expressions (->, <->, etc)
pub const InfixExpression = union(enum) {
    conditional: Conditional,
    biconditional: Biconditional,
    conjunction: Conjunction,
    disjunction: Disjunction,
};

pub const BlockExpression = struct {
    expressions: std.ArrayList(usize),
};

//an atom is a literal, something that fits in between connectives (e.g x,y,P)
pub const Atom = struct {
    symbol: Token,
};

//a predicate follows a quantifier (e.g P(x), P(x,y))
pub const Predicate = struct {
    //the variable literal right before the parens
    symbol: Token,
    //the list of args
    variables: std.ArrayList(Atom),
};

pub const Conjunction = struct {
    left: usize, //index to node
    right: usize, //index to node

    symbol: Token, //not sure if this is needed
};

pub const Disjunction = struct {
    left: usize, //index to node
    right: usize, //index to node

    symbol: Token, //not sure if this is needed
};

pub const Negation = struct {
    right: usize, //expression right of neg symbol

    symbol: Token, //might not be needed?
};

pub const Biconditional = struct {
    left: usize, //expression left of <->
    right: usize, //expression right of <->

    symbol: Token, //might not be needed?
};

pub const Conditional = struct {
    left: usize, //expression left of ->
    right: usize, //expression right of ->

    symbol: Token, //might not be needed?
};
