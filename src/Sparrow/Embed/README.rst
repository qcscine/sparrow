============================
SCINE - Sparrow Embed Binary
============================

What?
=====

This folder exists to embed Turbomole basis files and semi-empirical parameter
sets for NDDO and DFTB into C++ code.

Why?
====

It is necessary for these parameters to be hard-embedded into the binary in
order to ease binary relocatability of Sparrow. It is otherwise difficult to
determine the relative path of any packaged parameter files at runtime
correctly in all circumstances (loaded by app, loaded by module, loaded by
Python).

How? (practical)
================

Embed recognizes different use cases directly by its arguments.

- Directory: Assumed to contain a DFTB SKF parameter set
- ``*.basis``: Assumed to be a Turbomole format basis file
- ``*.json``: Assumed to be our in-house nddo parameter archive format

Embed generates the relevant code and places it next to its arguments.

You can pass as many arguments as you like.

How? (technical)
================

The binary itself is essentially a meta-compiler. It transforms its arguments
and sub-parameters into valid function identifiers and generates functions
generating runtime instances. It is not a ``constexpr`` meta-compiler because the
datatypes for the parameter sets are principally of variable length and
therefore do not lend themselves to ``constexpr`` embedding without addition of
runtime datatype conversion functions, significantly inflating binary size.

Possible improvements
=====================

As it stands, the ``embed`` binary itself depends on Sparrow. That's a little
awkward, since if we wanted a clean meta-compilation pipeline, it should be
upstream of Sparrow. The way it is, you need to compile Sparrow, then Embed,
then generate new parameter sets for sparrow, and recompile Sparrow.

To make a clean pipeline, separate out the runtime parameter types and their
parsers into a separate library and make Sparrow and Embed each depend on it.
Make sure Embed doesn't depend on Sparrow. Then you could integrate
meta-compilation with Embed into the build of Sparrow.

The way things stand are additionally impractical in that developers may easily
forget to update the code generating the runtime values for relocatability by
running the Embed binary again and shuffling the files into their proper
places. This would also be fixed by properly adding the meta-compilation step
into the build.
