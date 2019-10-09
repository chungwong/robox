defmodule Robox do
  alias Robox.Gbp1d

  @moduledoc """
  """

  @data %{
    oid: [1_428_571, 1_428_571, 1_428_571, 1_428_572, 1_428_572, 1_428_572, 1_428_572, 1_428_572],
    sku: ["A0A0A0", "A0A0A1", "A0A0A1", "A0A0A0", "A0A0A1", "A0A0A1", "A0A0A2", "A0A0A3"],
    l: [2.140000, 7.240000, 7.240000, 2.140000, 7.240000, 7.240000, 6.000000, 4.000000],
    d: [3.580000, 7.240000, 7.240000, 3.580000, 7.240000, 7.240000, 6.000000, 4.000000],
    h: [4.760000, 2.580000, 2.580000, 4.760000, 2.580000, 2.580000, 6.000000, 4.000000],
    # w:   [243.0000, 110.0000, 110.0000, 243.0000, 110.0000, 110.0000, 235.0000, 258.0000]
    w: [2.0000, 1.0000, 1.0000, 2.0000, 1.0000, 1.0000, 2.0000, 2.0000]
  }

  @doc """
  p: vector of items profit
  w: vector of items weight
  c: weight limit constraint
  k: indicator of which items are selected
  o: weight of items selected
  ok: indicator of whether all items are selected
  """
  @spec packager() :: nil
  def packager() do
    # mat = Matrex.new([@data.l |> Enum.slice(-5..-1), @data.d |> Enum.slice(-5..-1), @data.h |> Enum.slice(-5..-1)])
    p = Matrex.new([[36.46731, 135.23741, 135.23741, 216.00000, 64.00000]])
    w = Matrex.new([[2.0000, 1.0000, 1.0000, 2.0000, 2.0000]])

    gbp1d_solver_dpp(p, w, 15.28)
  end

  @doc """
  Return a column vector containing the indices of elements of X that are non-zero or satisfy a relational condition
  http://arma.sourceforge.net/docs.html#find
  """
  def find(mat, val) do
    Matrex.to_list(mat)
    |> Enum.with_index(1)
    |> Enum.flat_map(fn {m, i} ->
      if m == val do
        [i]
      else
        []
      end
    end)
  end

  @doc """
  For vector V, return true if all elements of the vector are non-zero or satisfy a relational condition
  http://arma.sourceforge.net/docs.html#all
  """
  @spec all?(struct, number) :: boolean
  def all?(mat, val) do
    Matrex.to_list(mat)
    |> Enum.all?(fn m ->
      m == val
    end)
  end

  def gbp1d_solver_dpp(p, w, c) do
    c = trunc(c)

    {_, n} = w[:size]

    k = Matrex.zeros(n, 1)

    q = Matrex.zeros(c + 1, n)
    g = Matrex.zeros(c + 1, 1)

    # forward
    {q, g} =
      Enum.reduce(1..n, {q, g}, fn i, {q, g} ->
        q = Matrex.set_column(q, i, g)

        w_i = trunc(w[i])

        g =
          Enum.reduce((w_i - 1)..(c - 1), g, fn j, g ->
            row = c - j + w_i

            val = max(g[row], g[c - j] + p[i])
            Matrex.set(g, row, 1, val)
          end)

        {q, g}
      end)

    qmax = g[c]

    # backward
    j = c

    {k, _, _} =
      Enum.reduce(1..n, {k, j, qmax}, fn i, {k, j, qptr} ->
        if q[j][n + 1 - i] < qptr do
          k = Matrex.set(k, n + 1 - i, 1, 1)
          j = trunc(j - w[n + 1 - i])
          qptr = q[j][n + 1 - i]

          {k, j, qptr}
        else
          {k, j, qptr}
        end
      end)

    s = find(k, 1)

    o =
      p[List.first(s)..List.last(s)]
      |> Enum.sum()

    ok = all?(k, 1)

    %Gbp1d{
      p: p,
      w: w,
      c: c,
      k: k,
      o: o,
      ok: ok
    }
  end
end
