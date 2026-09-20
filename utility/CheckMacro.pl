#!/usr/bin/perl
# -----------------------------------------------------------------------------
# CheckMacro.pl
#
# Add explicit compile-time checks for macros that are used as values in C/C++
# preprocessor expressions.  The C and C++ preprocessors normally replace an
# undefined identifier in a #if expression by zero.  That behavior can silently
# select an unintended AMPS configuration branch.  This utility converts that
# silent fallback into a clear #error diagnostic.
#
# AMPS invokes the legacy interface from its top-level Makefile:
#
#   ./utility/CheckMacro.pl build -in-place
#
# The interface above is intentionally preserved.  The implementation below is
# not a C++ parser, but it does implement the parts of preprocessing needed for
# this task:
#
#   * backslash-continued logical directives;
#   * stateful removal of // and /* */ comments;
#   * lexical recognition of identifiers, numbers, strings and characters;
#   * #if/#elif/#else/#endif nesting;
#   * correct association of #elif macros with their opening conditional;
#   * idempotent, marked generated blocks;
#   * atomic file replacement.
#
# The script is designed to operate on AMPS's disposable build tree.  It can be
# used on another tree, but --dry-run or --check-only should be used first.
# -----------------------------------------------------------------------------

use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Spec;
use File::Temp qw(tempfile);
use Getopt::Long qw(GetOptionsFromArray Configure);

my $BEGIN_MARKER = '/* AMPS_CHECKMACRO_BEGIN */';
my $END_MARKER   = '/* AMPS_CHECKMACRO_END */';

my %DEFAULT_EXTENSION = map { $_ => 1 }
    qw(h hpp c cc cpp cxx cu dfn);

my %DEFAULT_EXCLUDED_DIRECTORY = map { $_ => 1 }
    qw(.git .svn CVS);

# C++ alternative operator spellings are lexical operators, not macro names.
# true/false/nullptr are also excluded so that valid C++ preprocessor constants
# are not turned into #ifndef checks.
my %NON_MACRO_IDENTIFIER = map { $_ => 1 } qw(
    and and_eq bitand bitor compl false not not_eq nullptr or or_eq true xor
    xor_eq
);

# These compiler-provided preprocessor query operators accept tokens that are
# not necessarily macros (for example, an attribute name).  Their complete
# parenthesized argument list is therefore excluded from macro extraction.
my %PREPROCESSOR_QUERY_OPERATOR = map { $_ => 1 } qw(
    __has_attribute __has_builtin __has_cpp_attribute __has_extension
    __has_feature __has_include __has_include_next __has_warning
);

my %Option = (
    in_place  => 0,
    check_only => 0,
    dry_run   => 0,
    verbose   => 0,
);
my @ExcludedDirectory;
my $ExtensionList;
my $HeaderFile;
my $Help;

# Getopt::Long accepts the conventional --in-place form.  Normalize the two
# historical single-dash AMPS spellings before option processing.
for my $argument (@ARGV) {
    $argument = '--in-place' if $argument eq '-in-place';
    $argument = '--help'     if $argument eq '-help';
}

Configure(qw(no_auto_abbrev no_ignore_case permute));
GetOptionsFromArray(
    \@ARGV,
    'in-place!'    => \$Option{in_place},
    'check-only!'  => \$Option{check_only},
    'dry-run!'     => \$Option{dry_run},
    'header=s'     => \$HeaderFile,
    'verbose+'     => \$Option{verbose},
    'extensions=s' => \$ExtensionList,
    'exclude=s@'   => \@ExcludedDirectory,
    'help|h'       => \$Help,
) or usage(2);

usage(0) if $Help;
usage(2, 'Exactly one input directory is required') unless @ARGV == 1;

my $InputDirectory = abs_path($ARGV[0]);
die "CheckMacro.pl: input directory '$ARGV[0]' does not exist\n"
    unless defined($InputDirectory) && -d $InputDirectory;

die "CheckMacro.pl: --in-place and --header cannot be used together\n"
    if $Option{in_place} && defined($HeaderFile);

# With no explicit action, retain the legacy behavior of writing MacroCheck.h.
# Diagnostic modes analyze the in-place transformation without changing files.
if (!$Option{in_place} && !$Option{check_only} && !$Option{dry_run}
        && !defined($HeaderFile)) {
    $HeaderFile = 'MacroCheck.h';
}

my %Extension = %DEFAULT_EXTENSION;
if (defined($ExtensionList)) {
    %Extension = ();
    for my $extension (split(/,/, $ExtensionList)) {
        $extension =~ s/^\.//;
        die "CheckMacro.pl: empty extension in --extensions\n"
            if $extension eq '';
        $Extension{lc($extension)} = 1;
    }
}

my %ExcludedDirectory = %DEFAULT_EXCLUDED_DIRECTORY;
for my $directory (@ExcludedDirectory) {
    $ExcludedDirectory{$directory} = 1;
}

my @SourceFiles;
collect_source_files($InputDirectory, \@SourceFiles);

my %GlobalMacroSeen;
my @GlobalMacros;
my %Statistics = (
    files_scanned       => 0,
    files_with_checks   => 0,
    files_changed       => 0,
    conditional_groups  => 0,
    macro_uses          => 0,
);

for my $file (@SourceFiles) {
    my $analysis = analyze_file($file);
    ++$Statistics{files_scanned};
    $Statistics{conditional_groups} += scalar(@{$analysis->{groups}});
    $Statistics{macro_uses} += scalar(@{$analysis->{required_macros}});
    ++$Statistics{files_with_checks} if @{$analysis->{required_macros}};

    for my $macro (@{$analysis->{required_macros}}) {
        next if $GlobalMacroSeen{$macro}++;
        push(@GlobalMacros, $macro);
    }

    my $would_change = $analysis->{original_text} ne $analysis->{output_text};
    next unless $would_change;

    ++$Statistics{files_changed};
    print "CheckMacro.pl: would update $file\n"
        if $Option{dry_run} || $Option{check_only} || $Option{verbose};

    if ($Option{in_place} && !$Option{dry_run} && !$Option{check_only}) {
        atomic_write($file, $analysis->{output_text}, $analysis->{mode});
    }
}

if (defined($HeaderFile)) {
    my $header = generate_header(\@GlobalMacros);
    my $old_header = -e $HeaderFile ? read_file($HeaderFile) : '';

    if ($old_header ne $header) {
        print "CheckMacro.pl: would update $HeaderFile\n"
            if $Option{dry_run} || $Option{check_only} || $Option{verbose};

        if (!$Option{dry_run} && !$Option{check_only}) {
            my $mode = -e $HeaderFile ? (stat($HeaderFile))[2] : 0644;
            atomic_write($HeaderFile, $header, $mode);
        }
    }
}

if ($Option{verbose} || $Option{dry_run} || $Option{check_only}) {
    print "CheckMacro.pl: scanned $Statistics{files_scanned} files; " .
          "$Statistics{files_with_checks} contain checked expressions; " .
          "$Statistics{files_changed} would change; " .
          scalar(@GlobalMacros) . " unique macros found\n";
}

exit 0;

# =============================================================================
# Command-line help
# =============================================================================

sub usage {
    my ($exit_code, $message) = @_;
    print STDERR "CheckMacro.pl: $message\n\n" if defined($message);

    my $stream = $exit_code ? *STDERR : *STDOUT;
    print $stream <<'USAGE';
CheckMacro.pl adds explicit checks for macros used as values in C/C++
preprocessor conditional expressions.

Legacy AMPS usage:
  ./utility/CheckMacro.pl FOLDER -in-place
  ./utility/CheckMacro.pl FOLDER

Extended usage:
  ./utility/CheckMacro.pl [OPTIONS] FOLDER

Options:
  --in-place          Rewrite eligible source files atomically.
  --check-only        Analyze and report files that would change.
  --dry-run           Show the in-place changes without writing them.
  --header FILE       Write aggregate checks to FILE.  With no action option,
                      the legacy default is MacroCheck.h.
  --extensions LIST   Comma-separated extensions to scan.
  --exclude NAME      Exclude a directory basename; may be repeated.
  --verbose           Print processed-file and summary information.
  --help              Print this help text.

The default extensions are: h,hpp,c,cc,cpp,cxx,cu,dfn.
USAGE
    exit $exit_code;
}

# =============================================================================
# Deterministic directory traversal
# =============================================================================

sub collect_source_files {
    my ($directory, $files) = @_;

    opendir(my $handle, $directory)
        or die "CheckMacro.pl: cannot open directory '$directory': $!\n";
    my @entries = sort grep { $_ ne '.' && $_ ne '..' } readdir($handle);
    closedir($handle)
        or die "CheckMacro.pl: cannot close directory '$directory': $!\n";

    for my $entry (@entries) {
        my $path = File::Spec->catfile($directory, $entry);

        # Never follow directory symlinks: a recursive link can otherwise make
        # traversal unbounded or cause files outside the requested tree to be
        # rewritten.
        if (-d $path && !-l $path) {
            next if $ExcludedDirectory{$entry};
            collect_source_files($path, $files);
            next;
        }

        next unless -f $path;
        next if $entry =~ /~$/ || $entry =~ /^#.*#$/;
        next unless $entry =~ /\.([^.]+)$/;
        next unless $Extension{lc($1)};
        push(@$files, $path);
    }
}

# =============================================================================
# File analysis and conditional-group tracking
# =============================================================================

sub analyze_file {
    my ($file) = @_;
    my $mode = (stat($file))[2];
    my $original_text = read_file($file);

    # Generated blocks are removed before parsing.  This makes the operation
    # idempotent and prevents the checker's own #ifndef directives from being
    # interpreted as application conditionals on the next run.
    my $base_text = remove_generated_blocks($original_text, $file);
    my $clean_text = sanitize_comments($base_text);
    my @source_lines = split_lines($base_text);
    my @clean_lines  = split_lines($clean_text);

    die "CheckMacro.pl: internal line-map mismatch while processing '$file'\n"
        unless @source_lines == @clean_lines;

    my @stack;
    my @groups;

    for (my $line_index = 0; $line_index < @clean_lines; ++$line_index) {
        # Logical-line assembly is applied before recognizing a directive.  It
        # also skips continuation lines belonging to a multiline #define so a
        # '#' appearing there cannot be mistaken for a new directive.
        my $first_line = $line_index;
        my ($logical_line, $last_line) =
            assemble_logical_line(\@clean_lines, $first_line);
        $line_index = $last_line;

        next unless $logical_line =~ /^\s*#\s*([A-Za-z_][A-Za-z0-9_]*)\b(.*)$/s;
        my $directive = lc($1);
        my $argument = $2;

        if ($directive eq 'if' || $directive eq 'ifdef'
                || $directive eq 'ifndef') {
            my $group = {
                start_index => $first_line,
                macros => [],
                seen => {},
            };

            # The expression form requires lexical analysis.  #ifdef and
            # #ifndef intentionally ask whether a macro exists and therefore
            # must not be converted into mandatory-definition checks.
            if ($directive eq 'if') {
                add_group_macros($group, extract_required_macros($argument));
            }

            push(@stack, $group);
            next;
        }

        if ($directive eq 'elif') {
            die "$file:" . ($line_index + 1) .
                ": #elif without a matching opening conditional\n"
                unless @stack;
            add_group_macros($stack[-1], extract_required_macros($argument));
            next;
        }

        if ($directive eq 'else') {
            die "$file:" . ($line_index + 1) .
                ": #else without a matching opening conditional\n"
                unless @stack;
            next;
        }

        if ($directive eq 'endif') {
            die "$file:" . ($line_index + 1) .
                ": #endif without a matching opening conditional\n"
                unless @stack;
            push(@groups, pop(@stack));
            next;
        }
    }

    if (@stack) {
        my $line_number = $stack[-1]{start_index} + 1;
        die "$file:$line_number: opening conditional has no matching #endif\n";
    }

    # Groups close from the inside out, so restore source order before
    # generating insertion blocks and reporting macros.
    @groups = sort { $a->{start_index} <=> $b->{start_index} } @groups;

    my %insertion;
    my %required_seen;
    my @required_macros;
    for my $group (@groups) {
        next unless @{$group->{macros}};
        push(@{$insertion{$group->{start_index}}}, $group);
        for my $macro (@{$group->{macros}}) {
            next if $required_seen{$macro}++;
            push(@required_macros, $macro);
        }
    }

    my @output;
    for (my $i = 0; $i < @source_lines; ++$i) {
        if (exists($insertion{$i})) {
            for my $group (@{$insertion{$i}}) {
                push(@output, generate_check_block(
                    $source_lines[$i], $group->{macros}));
            }
        }
        push(@output, $source_lines[$i]);
    }

    return {
        original_text   => $original_text,
        output_text     => join('', @output),
        required_macros => \@required_macros,
        groups          => \@groups,
        mode            => $mode,
    };
}

sub add_group_macros {
    my ($group, @macros) = @_;
    for my $macro (@macros) {
        next if $group->{seen}{$macro}++;
        push(@{$group->{macros}}, $macro);
    }
}

# =============================================================================
# Physical-to-logical line assembly and comment sanitization
# =============================================================================

sub split_lines {
    my ($text) = @_;
    return () if $text eq '';

    my @lines = split(/(?<=\n)/, $text, -1);
    pop(@lines) if @lines && $lines[-1] eq '';
    return @lines;
}

sub assemble_logical_line {
    my ($lines, $first_index) = @_;
    my $last_index = $first_index;
    my $logical = $lines->[$last_index];

    while ($logical =~ /\\[ \t]*(?:\r?\n)?\z/
            && $last_index + 1 < @$lines) {
        $logical =~ s/\\[ \t]*(?:\r?\n)?\z/ /;
        ++$last_index;
        $logical .= $lines->[$last_index];
    }

    return ($logical, $last_index);
}

# Replace comments with spaces while preserving every newline and every byte
# position.  Preserving positions makes the sanitized and original physical
# line arrays directly interchangeable.  Quote states prevent URLs or comment
# delimiters inside ordinary string/character literals from being stripped.
sub sanitize_comments {
    my ($text) = @_;
    my $output = '';
    my $state = 'normal';
    my $escaped = 0;

    for (my $i = 0; $i < length($text); ++$i) {
        my $character = substr($text, $i, 1);
        my $next = $i + 1 < length($text) ? substr($text, $i + 1, 1) : '';

        if ($state eq 'line_comment') {
            if ($character eq "\n") {
                $output .= $character;
                $state = 'normal';
            }
            else {
                $output .= ' ';
            }
            next;
        }

        if ($state eq 'block_comment') {
            if ($character eq '*' && $next eq '/') {
                $output .= '  ';
                ++$i;
                $state = 'normal';
            }
            elsif ($character eq "\n" || $character eq "\r") {
                $output .= $character;
            }
            else {
                $output .= ' ';
            }
            next;
        }

        if ($state eq 'single_quote' || $state eq 'double_quote') {
            $output .= $character;
            if ($escaped) {
                $escaped = 0;
            }
            elsif ($character eq '\\') {
                $escaped = 1;
            }
            elsif (($state eq 'single_quote' && $character eq "'")
                    || ($state eq 'double_quote' && $character eq '"')) {
                $state = 'normal';
            }
            next;
        }

        if ($character eq '/' && $next eq '/') {
            $output .= '  ';
            ++$i;
            $state = 'line_comment';
        }
        elsif ($character eq '/' && $next eq '*') {
            $output .= '  ';
            ++$i;
            $state = 'block_comment';
        }
        elsif ($character eq "'") {
            $output .= $character;
            $state = 'single_quote';
            $escaped = 0;
        }
        elsif ($character eq '"') {
            $output .= $character;
            $state = 'double_quote';
            $escaped = 0;
        }
        else {
            $output .= $character;
        }
    }

    return $output;
}

# =============================================================================
# Preprocessor-expression lexer and macro extraction
# =============================================================================

sub extract_required_macros {
    my ($expression) = @_;
    my @tokens = lex_expression($expression);
    my %skip_token;
    my %defined_operand;

    # First identify `defined NAME` and `defined(NAME)`.  A macro explicitly
    # tested with defined is optional by construction.  All occurrences of that
    # name in the same expression are excluded, which supports the common safe
    # idiom `defined(FEATURE) && FEATURE == 1`.
    for (my $i = 0; $i < @tokens; ++$i) {
        next unless $tokens[$i]{type} eq 'identifier';
        next unless $tokens[$i]{text} eq 'defined';
        $skip_token{$i} = 1;

        my $j = $i + 1;
        if ($j < @tokens && $tokens[$j]{text} eq '(') {
            $skip_token{$j} = 1;
            ++$j;
        }
        if ($j < @tokens && $tokens[$j]{type} eq 'identifier') {
            $defined_operand{$tokens[$j]{text}} = 1;
            $skip_token{$j} = 1;
        }
    }

    # Compiler query operators are not ordinary function-like macros.  Skip the
    # operator and its complete argument list so attribute/header/builtin names
    # do not become false AMPS configuration requirements.
    for (my $i = 0; $i < @tokens; ++$i) {
        next unless $tokens[$i]{type} eq 'identifier';
        next unless $PREPROCESSOR_QUERY_OPERATOR{$tokens[$i]{text}};
        $skip_token{$i} = 1;
        next unless $i + 1 < @tokens && $tokens[$i + 1]{text} eq '(';

        my $depth = 0;
        for (my $j = $i + 1; $j < @tokens; ++$j) {
            $skip_token{$j} = 1;
            ++$depth if $tokens[$j]{text} eq '(';
            --$depth if $tokens[$j]{text} eq ')';
            if ($depth == 0) {
                $i = $j;
                last;
            }
        }
    }

    my %seen;
    my @macros;
    for (my $i = 0; $i < @tokens; ++$i) {
        next if $skip_token{$i};
        next unless $tokens[$i]{type} eq 'identifier';
        my $identifier = $tokens[$i]{text};
        next if $NON_MACRO_IDENTIFIER{$identifier};
        next if $defined_operand{$identifier};
        next if $seen{$identifier}++;
        push(@macros, $identifier);
    }

    return @macros;
}

sub lex_expression {
    my ($expression) = @_;
    my @tokens;
    my $length = length($expression);

    for (my $i = 0; $i < $length;) {
        my $character = substr($expression, $i, 1);

        if ($character =~ /\s/) {
            ++$i;
            next;
        }

        if ($character =~ /[A-Za-z_]/) {
            my $start = $i++;
            ++$i while $i < $length
                && substr($expression, $i, 1) =~ /[A-Za-z0-9_]/;
            push(@tokens, {
                type => 'identifier',
                text => substr($expression, $start, $i - $start),
            });
            next;
        }

        # Consume a complete preprocessing number before looking for
        # identifiers.  This prevents the xFF in 0xFF and suffixes such as UL
        # from being reported as macro names.
        if ($character =~ /[0-9]/
                || ($character eq '.' && $i + 1 < $length
                    && substr($expression, $i + 1, 1) =~ /[0-9]/)) {
            my $start = $i++;
            while ($i < $length) {
                my $current = substr($expression, $i, 1);
                my $previous = substr($expression, $i - 1, 1);
                if ($current =~ /[A-Za-z0-9_\.']/
                        || (($current eq '+' || $current eq '-')
                            && $previous =~ /[eEpP]/)) {
                    ++$i;
                }
                else {
                    last;
                }
            }
            push(@tokens, {
                type => 'number',
                text => substr($expression, $start, $i - $start),
            });
            next;
        }

        if ($character eq "'" || $character eq '"') {
            my $quote = $character;
            my $start = $i++;
            my $escaped = 0;
            while ($i < $length) {
                my $current = substr($expression, $i++, 1);
                if ($escaped) {
                    $escaped = 0;
                }
                elsif ($current eq '\\') {
                    $escaped = 1;
                }
                elsif ($current eq $quote) {
                    last;
                }
            }
            push(@tokens, {
                type => $quote eq "'" ? 'character' : 'string',
                text => substr($expression, $start, $i - $start),
            });
            next;
        }

        # Operators and punctuation do not need semantic classification for
        # macro extraction.  Retaining each character is sufficient for
        # recognizing parentheses around defined/query-operator operands.
        push(@tokens, { type => 'punctuation', text => $character });
        ++$i;
    }

    return @tokens;
}

# =============================================================================
# Generated output
# =============================================================================

sub generate_check_block {
    my ($opening_line, $macros) = @_;
    my ($indent) = $opening_line =~ /^(\s*)/;
    $indent = '' unless defined($indent);
    $indent =~ s/[\r\n]//g;
    my $newline = $opening_line =~ /\r\n/ ? "\r\n" : "\n";

    my $output = $indent . $BEGIN_MARKER . $newline;
    for my $macro (@$macros) {
        $output .= $indent . "#ifndef $macro" . $newline;
        $output .= $indent .
            "#error \"AMPS: $macro is used in a preprocessor expression " .
            "but is not defined\"" . $newline;
        $output .= $indent . '#endif' . $newline;
    }
    $output .= $indent . $END_MARKER . $newline;
    return $output;
}

sub remove_generated_blocks {
    my ($text, $file) = @_;
    my @lines = split_lines($text);
    my @output;
    my $inside = 0;

    for my $line (@lines) {
        if ($line =~ /^\s*\/\*\s*AMPS_CHECKMACRO_BEGIN\s*\*\/\s*(?:\r?\n)?\z/) {
            die "CheckMacro.pl: nested generated block in '$file'\n" if $inside;
            $inside = 1;
            next;
        }
        if ($line =~ /^\s*\/\*\s*AMPS_CHECKMACRO_END\s*\*\/\s*(?:\r?\n)?\z/) {
            die "CheckMacro.pl: unmatched generated-block end in '$file'\n"
                unless $inside;
            $inside = 0;
            next;
        }
        push(@output, $line) unless $inside;
    }

    die "CheckMacro.pl: unterminated generated block in '$file'\n" if $inside;
    return join('', @output);
}

sub generate_header {
    my ($macros) = @_;
    my $output = <<'HEADER';
/* Generated by utility/CheckMacro.pl.  Do not edit manually. */
#ifndef AMPS_GENERATED_MACRO_CHECK_H
#define AMPS_GENERATED_MACRO_CHECK_H

HEADER

    for my $macro (@$macros) {
        $output .= "#ifndef $macro\n";
        $output .= "#error \"AMPS: $macro is used in a preprocessor " .
                   "expression but is not defined\"\n";
        $output .= "#endif\n";
    }

    $output .= "\n#endif /* AMPS_GENERATED_MACRO_CHECK_H */\n";
    return $output;
}

# =============================================================================
# Checked and atomic file I/O
# =============================================================================

sub read_file {
    my ($file) = @_;
    open(my $handle, '<', $file)
        or die "CheckMacro.pl: cannot read '$file': $!\n";
    binmode($handle);
    local $/;
    my $content = <$handle>;
    close($handle)
        or die "CheckMacro.pl: cannot close '$file' after reading: $!\n";
    return defined($content) ? $content : '';
}

sub atomic_write {
    my ($file, $content, $mode) = @_;
    my $directory = dirname($file);
    my ($handle, $temporary) = tempfile(
        '.checkmacro.XXXXXX',
        DIR => $directory,
        UNLINK => 0,
    );

    eval {
        binmode($handle);
        print {$handle} $content
            or die "cannot write temporary file '$temporary': $!";
        close($handle)
            or die "cannot close temporary file '$temporary': $!";
        chmod($mode & 07777, $temporary)
            or die "cannot preserve permissions on '$temporary': $!";
        rename($temporary, $file)
            or die "cannot replace '$file' with '$temporary': $!";
    };

    if ($@) {
        my $error = $@;
        close($handle) if defined(fileno($handle));
        unlink($temporary) if -e $temporary;
        die "CheckMacro.pl: $error\n";
    }
}
