This is a direct port of https://github.com/gyang274/gbp

Reference [Paper](https://arxiv.org/abs/1809.10210)

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
```


## Installation
```elixir
def deps do
  [
    {:robox, "~> 0.1.0"}
  ]
end
```
Please refer to [Matrex](https://github.com/versilov/matrex) for OS dependencies for compiling
