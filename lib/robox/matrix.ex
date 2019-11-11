defmodule Robox.Matrix do
  def log(m) do
    Matrex.apply(m, :log)
  end

  @spec sum(%Matrex{}) :: number
  def sum(%Matrex{} = m) do
    Matrex.sum(m)
  end

  @spec sum(any) :: number
  def sum(_) do
    0.0
  end

  @doc """
  Fills a matrex with val at indices
  """
  def fill(m, nil, _) do
    m
  end

  @spec fill(%Matrex{}, %Matrex{}, number) :: %Matrex{}
  def fill(m, vector, val) do
    Enum.reduce(1..m[:cols], m, fn j, m ->
      Enum.reduce(Matrex.to_list(vector), m, fn i, m ->
        Matrex.set(m, trunc(i), j, val)
      end)
    end)
  end

  def prod(m, _dim \\ 0) do
    if m[:rows] == 1 || m[:cols] == 1 do
      Matrex.to_list(m)
      |> Enum.reduce(1, fn i, sum ->
        sum * i
      end)
    else
      throw("prod only supports vectors")
    end
  end

  @doc """
  m is a vector, index is one-based
  """
  @spec subvec(%Matrex{}, integer, integer) :: %Matrex{}
  def subvec(m, first_index, last_index)
      when is_integer(first_index) and is_integer(last_index) do
    list =
      Matrex.to_list(m)
      |> Enum.slice((first_index - 1)..(last_index - 1))

    # determine cols and rows of vector
    if m[:rows] > m[:cols] do
      Enum.map(list, fn l ->
        [l]
      end)
    else
      [list]
    end
    |> Matrex.new()
  end

  def boolean_to_integer(true), do: 1
  def boolean_to_integer(false), do: 0

  def column_to_list(m) do
    Enum.flat_map(1..m[:cols], fn i ->
      Matrex.column_to_list(m, i)
    end)
  end

  def concat(nil, m2) do
    m2
  end

  def concat(m1, nil) do
    m1
  end

  def concat(m1, m2) do
    Matrex.concat(m1, m2)
  end

  def linspace(_, _, n \\ 0)

  def linspace(_, _, 0) do
    nil
  end

  def linspace(0, 0, 1) do
    Matrex.new([[0]])
  end

  def linspace(_, endVal, 1) do
    Matrex.new([[endVal]])
  end

  def linspace(startVal, endVal, n) do
    n = trunc(n)
    step = (endVal - startVal) / (n - 1)

    Enum.map(0..(n - 1), fn i ->
      [startVal + step * i]
    end)
    |> Matrex.new()
  end

  @doc """
  Return a column vector containing the indices of elements of X that are non-zero or satisfy a relational condition
  http://arma.sourceforge.net/docs.html#find
  """
  @spec find({%Matrex{}, atom, number}, atom, {%Matrex{}, atom, number}) :: %Matrex{}
  def find({mat1, op1, val1} = c1, operator, {mat2, op2, val2} = c2)
      when is_tuple(c1) and is_tuple(c2) do
    rst1 =
      find(mat1, op1, val1)
      |> Matrex.to_list()

    rst2 =
      find(mat2, op2, val2)
      |> Matrex.to_list()

    case operator do
      :and ->
        Enum.reduce(rst1, [], fn i, rst ->
          Enum.reduce(rst2, rst, fn j, rst ->
            if i != nil && j != nil && i == j do
              List.insert_at(rst, -1, [i])
            else
              rst
            end
          end)
        end)

      :or ->
        (rst1 ++ rst2)
        |> Enum.uniq()
        |> Enum.sort()
        |> Enum.map(fn i -> [i] end)
    end
    |> case do
      [] -> nil
      m -> Matrex.new(m)
    end
  end

  def find(mat, op, val) do
    rst =
      Matrex.to_list(mat)
      |> Enum.with_index(1)

    Enum.reduce(rst, [], fn {m, i}, rst ->
      if op.(m, val) do
        List.insert_at(rst, -1, [i])
      else
        rst
      end
    end)
    |> Matrex.new()
  end

  def find(mat, val) do
    column_to_list(mat)
    |> Enum.with_index(1)
    |> Enum.flat_map(fn {m, i} ->
      if m == val do
        [i]
      else
        []
      end
    end)
    |> case do
      nil ->
        nil

      [] ->
        nil

      list ->
        Matrex.new([list])
    end
  end

  @doc """
  For vector V, return true if all elements of the vector are non-zero or satisfy a relational condition
  http://arma.sourceforge.net/docs.html#all
  """
  @spec all?(struct, number) :: boolean
  def all?(mat, val) do
    column_to_list(mat)
    |> Enum.all?(fn m ->
      m == val
    end)
  end

  @spec get_rows(%Matrex{}, nil) :: nil
  def get_rows(_, nil) do
    nil
  end

  @spec get_rows(%Matrex{}, []) :: nil
  def get_rows(_, []) do
    nil
  end

  @spec get_rows(%Matrex{}, list) :: %Matrex{}
  def get_rows(m, %Matrex{} = index) do
    get_rows(m, Matrex.to_list(index))
  end

  @spec get_rows(%Matrex{}, list) :: %Matrex{}
  def get_rows(m, index) when is_list(index) do
    Enum.map(index, fn i ->
      Matrex.row(m, trunc(i))
      |> Matrex.to_list()
    end)
    |> Matrex.new()
  end

  @spec get_columns(%Matrex{}, nil) :: nil
  def get_columns(_, nil) do
    nil
  end

  @spec get_columns(%Matrex{}, list) :: nil
  def get_columns(_, index) when index == [] do
    # FIXME
    nil
  end

  @spec get_columns(%Matrex{}, list) :: %Matrex{}
  def get_columns(m, %Matrex{} = index) do
    get_columns(m, Matrex.to_list(index))
  end

  @spec get_columns(%Matrex{}, list) :: %Matrex{}
  def get_columns(m, index) when is_list(index) do
    Enum.map(index, fn i ->
      Matrex.column(m, trunc(i))
    end)
    |> Matrex.concat()
  end

  @spec sort_index(list, [:asc | :desc]) :: %Matrex{}
  def sort_index(p, dir \\ :asc)

  def sort_index(%Matrex{} = m, dir) do
    m =
      case m[:rows] == 1 && m[:cols] > 1 do
        true -> Matrex.transpose(m)
        _ -> m
      end

    sort_index(column_to_list(m), dir)
  end

  def sort_index(list, dir) do
    sort_fun =
      if dir == :asc do
        &<=/2
      else
        &>=/2
      end

    Enum.with_index(list)
    |> Enum.sort_by(fn {v, _} -> v end, sort_fun)
    |> Enum.map(fn {_, i} -> [i] end)
    |> Matrex.new()
  end

  def stable_sort_index(list, dir \\ :asc) do
    sort_index(list, dir)
  end

  @spec set_rows(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def set_rows(m, index, vec) do
    # convert vector to list
    list = Matrex.to_list_of_lists(vec)

    list_index = Matrex.to_list(index)

    {m, _} =
      Enum.reduce(list_index, {m, list}, fn i, {m, list} ->
        i = trunc(i)
        {row, list} = List.pop_at(list, 0)

        m =
          Enum.reduce(1..length(row), m, fn j, m ->
            Matrex.set(m, i, j, Enum.at(row, j - 1))
          end)

        {m, list}
      end)

    m
  end

  @spec set_columns(%Matrex{}, nil, %Matrex{}) :: %Matrex{}
  def set_columns(m, nil, _) do
    m
  end

  @spec set_columns(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def set_columns(m, index, vec) do
    # convert vector to list
    list_index =
      Matrex.to_list(index)
      |> Enum.with_index(1)

    Enum.reduce(list_index, m, fn {e, i}, m ->
      e = trunc(e)
      Matrex.set_column(m, e, Matrex.column(vec, i))
    end)
  end

  @doc """
  http://arma.sourceforge.net/docs.html#sort
  """
  def sort(m, dir \\ :asc, dim \\ 0) do
    if m[:cols] == 1 and m[:rows] == 1 do
      throw("to sort")
    else
      if dim == 0 do
        Enum.map(1..m[:cols], fn i ->
          sort_column(Matrex.column(m, i), dir)
        end)
        |> Matrex.concat()
      else
        Enum.map(1..m[:rows], fn i ->
          sort_row(Matrex.row(m, i), dir)
        end)
        |> Matrex.new()
      end
    end
  end

  defp sort_column(%Matrex{} = col, dir) do
    [sort_column(Matrex.to_list(col), dir)]
    |> Matrex.new()
    |> Matrex.transpose()
  end

  defp sort_column(col, :desc) when is_list(col) do
    Enum.sort(col, &(&1 >= &2))
  end

  defp sort_column(col, :asc) when is_list(col) do
    Enum.sort(col)
  end

  defp sort_row(%Matrex{} = col, dir) do
    sort_column(Matrex.to_list(col), dir)
  end

  defp sort_row(row, :desc) when is_list(row) do
    Enum.sort(row, &(&1 >= &2))
  end

  defp sort_row(row, :asc) when is_list(row) do
    Enum.sort(row)
  end

  @doc """
  http://arma.sourceforge.net/docs.html#min_and_max
  """
  def max(m, dim \\ 0) do
    if dim == 0 do
      Enum.map(1..m[:cols], fn i ->
        [Matrex.max(Matrex.column(m, i))]
      end)
      |> Matrex.new()
      |> Matrex.transpose()
    else
      Enum.map(1..m[:rows], fn i ->
        [Matrex.max(Matrex.row(m, i))]
      end)
      |> Matrex.new()
    end
  end

  @doc """
  http://arma.sourceforge.net/docs.html#submat
  """
  @spec elem(%Matrex{}, %Matrex{}) :: %Matrex{}
  def elem(m, vector) do
    m = column_to_list(m)

    Matrex.to_list(vector)
    |> Enum.map(fn i ->
      [Enum.at(m, trunc(i - 1))]
    end)
    |> Matrex.new()
  end
end
