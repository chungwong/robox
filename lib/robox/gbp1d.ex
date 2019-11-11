defmodule Robox.Gbp1d do
  @moduledoc """
  p: vector of items profit
  w: vector of items weight
  c: weight limit constraint
  k: indicator of which items are selected
  o: weight of items selected
  ok: indicator of whether all items are selected
  """

  alias Robox.{Gbp1d, Matrix}
  defstruct [:p, :w, :c, :k, :o, :ok]

  def gbp1d_solver_dpp(p, w, c) do
    c = trunc(c)

    n = w[:cols]

    k = Matrex.zeros(n, 1)

    q = Matrex.zeros(c + 1, n)
    g = Matrex.zeros(c + 1, 1)

    # forward
    {q, g} =
      Enum.reduce(1..n, {q, g}, fn i, {q, g} ->
        q = Matrex.set_column(q, i, g)

        w_i = trunc(w[i])

        g =
          if w_i < c - 1 do
            Enum.reduce((w_i - 1)..(c - 1), g, fn j, g ->
              row = c - j + w_i

              val = max(Matrex.at(g, row, 1), Matrex.at(g, c - j, 1) + p[i])
              Matrex.set(g, row, 1, val)
            end)
          else
            g
          end

        {q, g}
      end)

    qmax = g[c + 1]

    # backward
    j = c + 1

    {k, _, _} =
      Enum.reduce(1..n, {k, j, qmax}, fn i, {k, j, qptr} ->
        if Matrex.at(q, j, n + 1 - i) < qptr do
          k = Matrex.set(k, n + 1 - i, 1, 1)
          j = trunc(j - w[n + 1 - i])
          qptr = Matrex.at(q, j, n + 1 - i)

          {k, j, qptr}
        else
          {k, j, qptr}
        end
      end)

    s = Matrix.find(k, 1)

    o =
      Matrix.get_columns(p, s)
      |> Matrix.sum()

    ok = Matrix.all?(k, 1)

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
