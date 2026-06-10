defmodule Robox do
  @moduledoc """
  Pure Elixir port of https://github.com/gyang274/gbp - generalized bin
  packing via extreme point heuristic and best information score fit.

      ldhw = [
        # [l, d, h, w]
        [2.14, 3.58, 4.76, 243.0],
        [7.24, 7.24, 2.58, 110.0],
        [6.0, 6.0, 6.0, 235.0],
        [4.0, 4.0, 4.0, 258.0]
      ]

      bn = [
        # [l, d, h, w]
        [6.0, 6.0, 6.0, 600.0]
      ]

      Robox.pack(ldhw, bn, :gbp4d)

  With a single bin, `:gbp4d` returns a `Robox.Gbp4d` struct; with several
  bins (sorted from smallest/most preferred to largest/least preferred) it
  returns a `Robox.Gbp4q` struct with the selected bin flagged in `f`.
  """

  alias Robox.{Gbp1d, Gbp4d, Gbp4q}

  @doc """
  Pack items into bin(s).

    * `:gbp4d` - `ldhw` items as `[l, d, h, w]` lists, `bn` a list of bins
      as `[l, d, h, w]` lists
    * `:gbp1d` - `ldhw` is `{p, w}` profits and weights, `bn` the capacity
  """
  @spec pack(list | {list, list}, list | number, atom) :: struct
  def pack(ldhw, bn, type \\ :gbp4d)

  def pack(ldhw, [bn], :gbp4d), do: Gbp4d.pack(ldhw, bn)
  def pack(ldhw, [_ | _] = bns, :gbp4d), do: Gbp4q.pack(ldhw, bns)
  def pack({p, w}, c, :gbp1d), do: Gbp1d.gbp1d_solver_dpp(p, w, c)
end
