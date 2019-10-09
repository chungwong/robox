defmodule Robox.Gbp4d.Ck do
  alias Robox.{Matrix}
  require Logger

  def gbp4d_checkr(sn) do
    klmt = Matrix.find(sn.k, 1)

    if klmt do
      row_or_col = if klmt[:rows] > klmt[:cols], do: :rows, else: :cols

      # main
      # main: no conflict between it and bn
      okfit =
        Enum.reduce(1..klmt[row_or_col], true, fn i, okfit ->
          ki = trunc(klmt[i])

          if Matrex.at(sn.it, 1, ki) + Matrex.at(sn.it, 5, ki) > sn.bn[1] ||
               Matrex.at(sn.it, 2, ki) + Matrex.at(sn.it, 6, ki) > sn.bn[2] ||
               Matrex.at(sn.it, 3, ki) + Matrex.at(sn.it, 7, ki) > sn.bn[3] do
            throw("gbp4d_checkr: it conflict bn: index #{ki}")
          else
            okfit
          end
        end)

      # main: no conflict between each pair of it
      okfit =
        Enum.reduce(1..klmt[row_or_col], okfit, fn i, okfit ->
          Enum.reduce((i + 1)..klmt[row_or_col], okfit, fn j, okfit ->
            if i < klmt[row_or_col] do
              ki = trunc(klmt[i])
              kj = trunc(klmt[j])

              if !(sn.it[1][ki] + sn.it[5][ki] <= sn.it[1][kj] ||
                     sn.it[1][kj] + sn.it[5][kj] <= sn.it[1][ki] ||
                     sn.it[2][ki] + sn.it[6][ki] <= sn.it[2][kj] ||
                     sn.it[2][kj] + sn.it[6][kj] <= sn.it[2][ki] ||
                     sn.it[3][ki] + sn.it[7][ki] <= sn.it[3][kj] ||
                     sn.it[3][kj] + sn.it[7][kj] <= sn.it[3][ki]) do
                throw("gbp4d_checkr: it conflict bn: index #{ki} and #{kj}")
              else
                okfit
              end
            else
              okfit
            end
          end)
        end)

      # main: no conflict on weight limit
      val =
        klmt
        |> Matrex.subtract(1)
        |> Matrex.multiply(8)
        |> Matrex.add(8)

      if Matrix.sum(Matrix.elem(sn.it, val)) > sn.bn[4] do
        throw("gbp4d_checkr: it conflict bn: conflict on weight constraint")
      else
        okfit
      end
    else
      true
    end
  end
end
