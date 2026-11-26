#!/usr/bin/env perl
# Observation score of the 2013-04-11 real-event srcSEP nightly test.
#
# Usage: score_observations.pl OBSERVATIONS.csv PLOTS_DIR OUTPUT_PREFIX MOVER
#
# Reads the comparison table produced at event setup (observations.csv) and
# the native AMPS observer products PLOTS_DIR/sample.flux/observer=ID.*.dat,
# writes OUTPUT_PREFIX.txt (summary) and OUTPUT_PREFIX.csv (row by row).
# The score is REPORTED, never gating: the nightly PASS/FAIL of the Table
# entry comes only from DiffNum.pl against the AMPS reference outputs.  The
# script therefore exits 0 after writing its report; a missing product is
# reported as status=ERROR inside the report rather than hidden.
#
# Model operator (srcSEP/sep.h cSamplingBuffer): a sample.flux column is the
# omnidirectional speed-weighted number density in one energy bin,
#   F = sum(W |v|) / V_segment                          [m^-2 s^-1 per bin],
# averaged over the output window [(k-1) C, k C) whose end k C is the time
# column (C = time_step * 12 output-cadence steps).  Hence
#   differential intensity  J  = F / (4 pi dE)          [(m^2 s sr J)^-1]
#                            -> x 1e-4 / (4 pi dE_MeV)   [(cm^2 s sr MeV)^-1]
#   integral intensity  J_int  = F / (4 pi) x 1e-4       [pfu]
#   fluence                    = 4 pi Int J_int dt = sum(F x 1e-4 x C) [cm^-2]
# (the last is the OpSEP fluence definition used by the CCMC GOES record).
#
# Calibration (user decision 2026-11-24): one global amplitude A from the
# rows with role=calibration only (SOHO/ERNE HED 4 h window),
#   ln A = mean( ln(obs) - ln(model) ),
# i.e. the least-squares amplitude in log space.  The source is linear in its
# injection efficiency when particle-wave coupling is off, so A rescales the
# run's nominal efficiency.  All role=holdout rows are scored with that A
# unchanged; role=excluded rows are listed but not scored.
# Metrics follow the existing OV01 definitions: point coverage (model and
# observation both positive), log10 RMSE, median |log10 error| and the
# Pearson correlation of log10 values, over all scored hold-out points.
use strict;
use warnings;
use POSIX qw(floor);
use List::Util qw(sum max);

my ($obs_csv, $plots, $prefix, $mover) = @ARGV;
die "usage: $0 OBSERVATIONS.csv PLOTS_DIR OUTPUT_PREFIX MOVER\n"
  unless defined $mover;
my $PI = 4.0 * atan2(1.0, 1.0);
my @report;
my $status = 'OK';

# --- observations ------------------------------------------------------------
open(my $in, '<', $obs_csv) or die "cannot read $obs_csv: $!\n";
my $header = <$in>;
chomp $header;
my @fields = split /,/, $header;
my @rows;
while (my $line = <$in>) {
  chomp $line;
  # The note column may contain commas only inside double quotes.
  my @values;
  while ($line =~ /\G("([^"]*(?:""[^"]*)*)"|[^,]*)(,|$)/g) {
    my $v = defined $2 ? $2 : $1;
    $v =~ s/""/"/g;
    push @values, $v;
    last if $3 eq '';
  }
  my %row;
  @row{@fields} = @values;
  push @rows, \%row;
}
close $in;

# --- model products ------------------------------------------------------------
my %series;  # observer id -> [ [t, F], ... ]
my $cadence;
sub load_observer {
  my ($id) = @_;
  return $series{$id} if exists $series{$id};
  my @files = glob("$plots/sample.flux/observer=$id.field-line=0.r=*.dat");
  if (@files != 1) {
    $series{$id} = undef;
    return undef;
  }
  open(my $f, '<', $files[0]) or return ($series{$id} = undef);
  my @data;
  while (my $line = <$f>) {
    next if $line =~ /^\s*(VARIABLES|ZONE|TITLE)/i;
    my @v = split ' ', $line;
    next unless @v >= 2;
    push @data, [ $v[0] + 0.0, $v[1] + 0.0 ];
  }
  close $f;
  if (@data >= 2 && !defined $cadence) { $cadence = $data[1][0] - $data[0][0]; }
  return ($series{$id} = \@data);
}

# Window mean of F over outputs whose interval [(t-C), t) lies in [a, b).
sub window_mean {
  my ($data, $a, $b) = @_;
  my @sel = grep { $_->[0] - $cadence >= $a - 1e-6 && $_->[0] <= $b + 1e-6 } @$data;
  return (undef, 0) unless @sel;
  return (sum(map { $_->[1] } @sel) / @sel, scalar @sel);
}

sub model_value {
  my ($row) = @_;
  my $data = load_observer($row->{observer});
  return (undef, 'missing observer product') unless defined $data && @$data >= 2;
  my ($a, $b) = ($row->{model_start_s}, $row->{model_stop_s});
  my $de = $row->{energy_high_mev} - $row->{energy_low_mev};
  if ($row->{quantity} eq 'differential') {
    my ($f, $n) = window_mean($data, $a, $b);
    return (undef, 'window outside model output') unless $n;
    return ($f * 1.0e-4 / (4.0 * $PI * $de), "$n outputs");
  }
  my @sel = grep { $_->[0] - $cadence >= $a - 1e-6 && $_->[0] <= $b + 1e-6 } @$data;
  return (undef, 'window outside model output') unless @sel;
  if ($row->{quantity} eq 'integral_peak') {
    return (max(map { $_->[1] } @sel) * 1.0e-4 / (4.0 * $PI), scalar(@sel) . ' outputs');
  }
  if ($row->{quantity} eq 'integral_fluence') {
    return (sum(map { $_->[1] } @sel) * 1.0e-4 * $cadence, scalar(@sel) . ' outputs');
  }
  return (undef, "unknown quantity $row->{quantity}");
}

for my $row (@rows) {
  my ($m, $note) = model_value($row);
  $row->{model} = $m;
  $row->{model_note} = $note;
}

# --- calibration -------------------------------------------------------------
my @calibration = grep { $_->{role} eq 'calibration' } @rows;
my @cal_ok = grep { defined $_->{model} && $_->{model} > 0 && $_->{value} > 0 } @calibration;
my $amplitude;
if (@cal_ok) {
  $amplitude = exp(sum(map { log($_->{value}) - log($_->{model}) } @cal_ok) / @cal_ok);
} else {
  $status = 'ERROR';
  push @report, 'ERROR: no calibration row has a positive model and observation';
}

# --- hold-out metrics ----------------------------------------------------------
sub metrics {
  my ($set) = @_;
  my @all = @$set;
  my @ok = grep { defined $_->{model} && $_->{model} > 0 && $_->{value} > 0 } @all;
  my %m = (points => scalar(@all), scored => scalar(@ok));
  return \%m unless @all;
  $m{coverage} = @ok / @all;
  return \%m unless @ok && defined $amplitude;
  my @err = map { (log($amplitude * $_->{model}) - log($_->{value})) / log(10.0) } @ok;
  $m{rmse_dex} = sqrt(sum(map { $_ * $_ } @err) / @err);
  my @abs = sort { $a <=> $b } map { abs $_ } @err;
  $m{median_abs_dex} = @abs % 2 ? $abs[$#abs / 2] : 0.5 * ($abs[@abs / 2 - 1] + $abs[@abs / 2]);
  $m{mean_bias_dex} = sum(@err) / @err;
  if (@ok >= 3) {
    my @x = map { log($amplitude * $_->{model}) / log(10.0) } @ok;
    my @y = map { log($_->{value}) / log(10.0) } @ok;
    my ($mx, $my) = (sum(@x) / @x, sum(@y) / @y);
    my ($sxy, $sxx, $syy) = (0, 0, 0);
    for my $i (0 .. $#x) {
      $sxy += ($x[$i] - $mx) * ($y[$i] - $my);
      $sxx += ($x[$i] - $mx) ** 2;
      $syy += ($y[$i] - $my) ** 2;
    }
    $m{correlation} = ($sxx > 0 && $syy > 0) ? $sxy / sqrt($sxx * $syy) : undef;
  }
  return \%m;
}

my @holdout = grep { $_->{role} eq 'holdout' } @rows;
my %groups = (
  'all hold-out' => \@holdout,
  'ERNE-HED 12 h' => [ grep { $_->{id} =~ /^erne-12h/ } @holdout ],
  'ERNE-HED 36 h' => [ grep { $_->{id} =~ /^erne-36h/ } @holdout ],
  'Liu Fig.12 12 h' => [ grep { $_->{id} =~ /^liu-12h/ } @holdout ],
  'Liu Fig.12 36 h' => [ grep { $_->{id} =~ /^liu-36h/ } @holdout ],
  'GOES-13 >10 MeV' => [ grep { $_->{id} =~ /^goes13/ } @holdout ],
  'calibration (in-sample)' => \@calibration,
);

sub fmt { my ($v, $f) = @_; return defined $v ? sprintf($f, $v) : 'n/a'; }
$status = 'ERROR' if grep { !defined $_->{model} } @holdout, @calibration;
push @report, "srcSEP 2013-04-11 real-event observation score (REPORTED, NOT GATING)";
push @report, "mover: $mover";
push @report, "status: $status";
push @report, "model output cadence [s]: " . fmt($cadence, '%.6g');
push @report, "calibration rows: " . scalar(@calibration) . " (usable " . scalar(@cal_ok) . ")";
push @report, "global amplitude A (multiplies nominal injection efficiency): " . fmt($amplitude, '%.6e');
for my $name (sort keys %groups) {
  my $m = metrics($groups{$name});
  push @report, sprintf("%-24s points=%3d scored=%3d coverage=%s rmse_dex=%s median_abs_dex=%s bias_dex=%s r=%s",
    $name, $m->{points}, $m->{scored}, fmt($m->{coverage}, '%.4f'), fmt($m->{rmse_dex}, '%.4f'),
    fmt($m->{median_abs_dex}, '%.4f'), fmt($m->{mean_bias_dex}, '%+.4f'), fmt($m->{correlation}, '%.4f'));
}
for my $row (grep { $_->{id} =~ /^goes13/ } @holdout) {
  push @report, sprintf("%s: observed=%.6g %s model(xA)=%s ratio=%s", $row->{id}, $row->{value}, $row->{unit},
    fmt(defined $row->{model} && defined $amplitude ? $amplitude * $row->{model} : undef, '%.6g'),
    fmt(defined $row->{model} && defined $amplitude && $row->{value} > 0 ? $amplitude * $row->{model} / $row->{value} : undef, '%.4f'));
}

open(my $txt, '>', "$prefix.txt") or die "cannot write $prefix.txt: $!\n";
print $txt "$_\n" for @report;
close $txt;
open(my $csv, '>', "$prefix.csv") or die "cannot write $prefix.csv: $!\n";
print $csv "id,role,instrument,quantity,energy_low_mev,energy_high_mev,model_start_s,model_stop_s,observed,unit,model_raw,model_scaled,log10_ratio,model_note\n";
for my $row (@rows) {
  my $scaled = (defined $row->{model} && defined $amplitude) ? $amplitude * $row->{model} : undef;
  my $ratio = (defined $scaled && $scaled > 0 && $row->{value} > 0) ? log($scaled / $row->{value}) / log(10.0) : undef;
  print $csv join(',', $row->{id}, $row->{role}, $row->{instrument}, $row->{quantity}, $row->{energy_low_mev},
    $row->{energy_high_mev}, $row->{model_start_s}, $row->{model_stop_s}, $row->{value}, $row->{unit},
    fmt($row->{model}, '%.6e'), fmt($scaled, '%.6e'), fmt($ratio, '%.4f'), $row->{model_note}), "\n";
}
close $csv;
print "$_\n" for @report;
exit 0;
