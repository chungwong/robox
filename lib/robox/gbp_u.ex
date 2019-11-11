defmodule Robox.GbpU do
  alias Robox.{Matrix}

  def create_max_mode_tbl(_m, _ulmt, _vlmt, _mval, _mlvl) do
  end

  @spec unique_rows(%Matrex{}) :: %Matrex{}
  def unique_rows(m) do
    ulmt =
      try do
        ulmt = Matrex.zeros(m[:rows], 1)

        Enum.each(1..m[:rows], fn i ->
          Enum.each((i + 1)..m[:rows], fn j ->
            if approx_equal_cpp(Matrex.row(m, i), Matrex.row(m, j)) do
              throw({:break, Matrex.set(ulmt, j, 1, 1)})
            end
          end)
        end)
      catch
        {:break, ulmt} ->
          ulmt
      end

    Matrix.get_columns(m, Matrix.find(ulmt, 0))
  end

  @spec unique_cols(%Matrex{}) :: %Matrex{}
  def unique_cols(m) do
    vlmt = Matrex.zeros(m[:cols], 1)

    vlmt =
      Enum.reduce(1..m[:cols], vlmt, fn i, vlmt ->
        start = if i + 1 > m[:cols], do: m[:cols], else: i + 1

        Enum.reduce_while(start..m[:cols], vlmt, fn j, vlmt ->
          if i != j && approx_equal_cpp(Matrex.column(m, i), Matrex.column(m, j)) do
            {:halt, Matrex.set(vlmt, j, 1, 1)}
          else
            {:cont, vlmt}
          end
        end)
      end)

    Matrix.get_columns(m, Matrix.find(vlmt, 0))
  end

  defp is_size_equal?(%Matrex{} = a, %Matrex{} = b) do
    Matrex.size(a) == Matrex.size(b)
  end

  @spec approx_equal_cpp(%Matrex{}, %Matrex{}, number) :: boolean
  def approx_equal_cpp(lhs, rhs, tol \\ 0.00000001) do
    approx_equal(lhs, rhs, :absdiff, tol)
  end

  @doc """
  Moves me to another module
  """
  def approx_equal(%Matrex{} = a, %Matrex{} = b, :absdiff, tol) do
    if is_size_equal?(a, b) do
      a_list = Matrix.column_to_list(a)
      b_list = Matrix.column_to_list(b)

      Enum.map(0..(length(a_list) - 1), fn i ->
        absdiff_approx_equal(Enum.at(a_list, i), Enum.at(b_list, i), tol)
      end)
      |> Enum.all?()
    else
      false
    end
  end

  def approx_equal(%Matrex{} = a, %Matrex{} = b, :reldiff, tol) do
    if is_size_equal?(a, b) do
      a_list = Matrix.column_to_list(a)
      b_list = Matrix.column_to_list(b)

      Enum.map(0..(length(a_list) - 1), fn i ->
        reldiff_approx_equal(Enum.at(a_list, i), Enum.at(b_list, i), tol)
      end)
      |> Enum.all?()
    else
      false
    end
  end

  def approx_equal(%Matrex{} = a, %Matrex{} = b, :both, tol) do
    approx_equal(a, b, :absdiff, tol) || approx_equal(a, b, :reldiff, tol)
  end

  def approx_equal(a, b, :both, tol) do
    approx_equal(a, b, :absdiff, tol) || approx_equal(a, b, :reldiff, tol)
  end

  def absdiff_approx_equal(a, b, tol) do
    if abs(a - b) <= tol do
      true
    else
      false
    end
  end

  def reldiff_approx_equal(a, b, tol) do
    if abs(a - b) / max(abs(a), abs(b)) <= tol do
      true
    else
      false
    end
  end

  def sort_index_via_cols_internal(_m, _ulmt, _vlmt) do
  end

  def sort_index_via_cols(_m, _vlmt) do
  end

  def sort_via_cols(_m, _vlmt) do
  end

  @spec sort_index_via_rows_internal(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def sort_index_via_rows_internal(m, ulmt, vlmt) do
    if m == nil || m[:rows] == 0 do
      []
    end

    try do
      if ulmt == nil do
        if vlmt == nil do
          throw({:return, nil})
          nil
        else
          throw({:return, Matrix.linspace(1, length(vlmt), length(vlmt))})
        end
      end

      mlmt =
        Matrex.multiply(vlmt, m[:rows])
        |> Matrex.add(ulmt[1])

      listm = Matrix.column_to_list(m)

      id =
        Enum.map(Matrix.column_to_list(mlmt), fn i ->
          Enum.at(listm, trunc(i))
        end)
        |> Matrix.stable_sort_index(:asc)

      if ulmt[:rows] == 1 do
        id
      else
        g = Matrex.fill(vlmt[:rows], 1, :nan)

        {_, g, _} =
          Enum.reduce(1..vlmt[:rows], {0, g, :nan}, fn i, {j, g, v} ->
            r = trunc(ulmt[1]) + 1
            c = trunc(vlmt[trunc(id[i]) + 1]) + 1

            if v == Matrex.at(m, r, c) do
              g = Matrex.set(g, i, 1, j)
              {j, g, v}
            else
              j = i
              g = Matrex.set(g, i, 1, j)
              v = Matrex.at(m, r, c)
              {j, g, v}
            end
          end)

        Enum.reduce(1..vlmt[:rows], id, fn i, id ->
          g0 = Matrix.find(g, i)

          if g0 && g0[:cols] > 1 do
            id0 = Matrix.get_rows(id, g0)

            id0id =
              sort_index_via_rows_internal(
                m,
                Matrix.subvec(ulmt, 2, ulmt[:rows]),
                Matrix.get_rows(vlmt, Matrex.add(id0, 1))
              )

            Matrix.set_rows(id, g0, Matrix.get_rows(id0, Matrex.add(id0id, 1)))
          else
            id
          end
        end)
      end
    catch
      {:return, vec} ->
        vec
    end
  end

  @spec sort_index_via_rows(%Matrex{}, %Matrex{}) :: %Matrex{}
  def sort_index_via_rows(m, ulmt) do
    cond do
      m == nil || m[:rows] == 0 || m[:cols] == 0 ->
        nil

      ulmt == nil ->
        Matrix.linspace(1, m[:cols], m[:cols])

      true ->
        sort_index_via_rows_internal(m, ulmt, Matrix.linspace(0, m[:cols] - 1, m[:cols]))
    end
  end

  @spec sort_via_rows(nil, %Matrex{}) :: nil
  def sort_via_rows(nil, _) do
    nil
  end

  @spec sort_via_rows(%Matrex{}, %Matrex{}) :: %Matrex{}
  def sort_via_rows(%Matrex{} = m, ulmt) do
    if m == nil || m[:rows] == 0 || m[:cols] == 0 || ulmt == nil do
      m
    else
      id =
        sort_index_via_rows_internal(m, ulmt, Matrix.linspace(0, m[:cols] - 1, m[:cols]))
        |> Matrex.add(1)

      Matrix.get_columns(m, id)
    end
  end
end
