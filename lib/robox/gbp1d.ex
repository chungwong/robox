defmodule Robox.Gbp1d do
  @moduledoc """
  Generalized bin packing problem in 1 dimension, a.k.a knapsack 0-1 problem.

  Solves

      maximize   sum p_j * k_j
      subject to sum w_j * k_j <= c, k_j in {0, 1}

  via dynamic programming, mirroring `gbp1d_solver_dpp` in gbp/src/gbp1d.cpp.

  p: vector of items profit
  w: vector of items weight (non-negative integers)
  c: weight limit constraint
  k: indicator of which items are selected
  o: profit of items selected
  ok: indicator of whether all items are selected
  """

  defstruct [:p, :w, :c, :k, :o, :ok]

  @type t :: %__MODULE__{
          p: [float],
          w: [non_neg_integer],
          c: non_neg_integer,
          k: [0 | 1],
          o: float,
          ok: boolean
        }

  @spec gbp1d_solver_dpp([number], [number], number) :: t
  def gbp1d_solver_dpp(p, w, c) do
    c = trunc(c)
    w = Enum.map(w, &trunc/1)
    p = Enum.map(p, &(&1 / 1))
    n = length(w)

    if n == 0 do
      %__MODULE__{p: p, w: w, c: c, k: [], o: 0.0, ok: true}
    else
      pt = List.to_tuple(p)
      wt = List.to_tuple(w)

      # forward: g is the DP value vector over capacities 0..c;
      # snapshots keep g before each item is processed, for backtracking
      g0 = List.duplicate(0.0, c + 1)

      {g, snapshots} =
        Enum.reduce(0..(n - 1), {g0, []}, fn i, {g, snaps} ->
          wi = elem(wt, i)
          pi = elem(pt, i)
          gt = List.to_tuple(g)

          g2 =
            Enum.map(0..c, fn cap ->
              v = elem(gt, cap)

              if cap >= wi do
                max(v, elem(gt, cap - wi) + pi)
              else
                v
              end
            end)

          {g2, [gt | snaps]}
        end)

      # backward: snapshots is [g_{n-1}, ..., g_0]
      qmax = List.last(g)

      {k, _, _} =
        Enum.reduce(Enum.with_index(snapshots), {[], c, qmax}, fn {snap, rev_i}, {ks, j, qptr} ->
          i = n - 1 - rev_i

          if elem(snap, j) < qptr do
            j2 = j - elem(wt, i)
            {[1 | ks], j2, elem(snap, j2)}
          else
            {[0 | ks], j, qptr}
          end
        end)

      o =
        Enum.zip(p, k)
        |> Enum.reduce(0.0, fn
          {pi, 1}, acc -> acc + pi
          {_, 0}, acc -> acc
        end)

      %__MODULE__{p: p, w: w, c: c, k: k, o: o, ok: Enum.all?(k, &(&1 == 1))}
    end
  end
end
