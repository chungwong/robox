# Robox

Pure Elixir port of [gbp](https://github.com/gyang274/gbp) - generalized bin
packing (3D bin packing with weight limit) via extreme point heuristic and
best information score fit strategy.

Reference [Paper](https://arxiv.org/abs/1809.10210)

No native dependencies - the previous Matrex-based implementation has been
rewritten in plain Elixir, which made it both orders of magnitude faster and
immune to the NIF crashes of the old version.

## Getting started

```elixir
ldhw =
  [
    # [l, d, h, w],
    [2.140000, 3.580000, 4.760000, 243.0000],
    [7.240000, 7.240000, 2.580000, 110.0000],
    [6.000000, 6.000000, 6.000000, 235.0000],
    [4.000000, 4.000000, 4.000000, 258.0000]
  ]

bn =
  [
    # [l, d, h, w],
    [06.0000, 06.0000, 06.0000, 600.000]
  ]

Robox.pack(ldhw, bn, :gbp4d)
#=> %Robox.Gbp4d{
#     it: [{x, y, z, w, l, d, h, wt}, ...],  # one per item, x = -1.0 when not fitted
#     k: [0 | 1, ...],                       # selection indicator per item
#     o: 550.47,                             # total fitted volume
#     ok: false,                             # did all items fit?
#     ...
#   }
```

With several bins (ordered from smallest / most preferred to largest / least
preferred), the most preferable bin that fits all (or most) items is selected
and flagged in `f`:

```elixir
bns = [
  [06.0, 06.0, 06.0, 600.0],
  [10.0, 08.0, 06.0, 600.0],
  [22.0, 14.0, 09.0, 800.0]
]

Robox.pack(ldhw, bns, :gbp4d)
#=> %Robox.Gbp4q{f: [0, 1, 0], ok: true, ...}
```

The 0-1 knapsack solver is also exposed:

```elixir
Robox.pack({[36.4, 135.2, 216.0], [243, 110, 235]}, 714, :gbp1d)
#=> %Robox.Gbp1d{k: [0, 1, 1], o: 351.2, ok: false}
```

A solution can be validated with `Robox.Gbp4d.Ck.gbp4d_checkr/1` /
`Robox.Gbp4d.Ck.gbp4q_checkr/1`.

## Correctness

`scripts/diff_vs_cpp.exs` runs randomized differential tests against the
original C++ implementation, comparing the profit ranking `p`, the selection
`k`, the objective `o` and `ok`. Placements can differ from the C++ output
when several placements score equally (armadillo's unstable sort breaks ties
arbitrarily; this port is deterministic), but the fitted set and objective
match.

## Installation

```elixir
def deps do
  [
    {:robox, "~> 0.2.0"}
  ]
end
```
