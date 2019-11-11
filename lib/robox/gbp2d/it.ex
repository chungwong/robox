defmodule Robox.Gbp2d.It do
  alias Robox.{Gbp2d.Ktlist2d}
  @spec gbp2d_it_create_ktlist(%Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}, number) :: struct
  def gbp2d_it_create_ktlist(bn, it, ktinit, xp, nlmt) do
    ktldht = gbp2d_it_create_ktldht(ktinit)
    kt = Matrex.zeros(4, ktldht[:cols] * xp[:cols])
    vlmt = Matrex.zeros(ktldht[:cols] * xp[:cols], 1)
    ij = 0
    {kt, vlmt} =
      Enum.reduce(1..ktldht[:cols], {kt, vlmt}, fn i, {kt, vlmt} ->
        Enum.reduce(1..xp[:cols], {kt, vlmt}, fn j, {kt, vlmt} ->
          if(Matrex.at(ktldht, 1, i) <= Matrex.at(xp, 3, j)
            && Matrex.at(ktldht, 2, i) <= Matrex.at(xp, 4, j)) do
            ij = (i - 1) * xp[:cols] + j # +1: Elixir index starts from 1 instead of 0
            kt = kt
                 |> Matrex.set(1, ij, Matrex.at(xp, 1, j))
                 |> Matrex.set(2, ij, Matrex.at(xp, 2, j))
                 |> Matrex.set(3, ij, Matrex.at(ktldht, 1, i))
                 |> Matrex.set(4, ij, Matrex.at(ktldht, 2, i))
            vlmt = Matrex.set(vlmt, ij, 1, 1)
            {kt, vlmt}
          else
            {kt, vlmt}
          end
        end)
      end)
    kt = Robox.Matrix.get_columns(kt, Robox.Matrix.find(vlmt, 1))
    kt_n_cols =
      if kt == nil do
        0
      else
        kt[:cols]
      end
    it_n_cols =
      if it == nil do
        0
      else
        it[:cols]
      end

    vlmt =
      if kt_n_cols != 0 do
        vlmt = Matrex.ones(kt_n_cols, 1)
        Enum.reduce(1..kt_n_cols, vlmt, fn i, vlmt ->
          if(it_n_cols != 0) do
            Enum.reduce_while(1..it_n_cols, vlmt, fn j, vlmt ->
              if !(Matrex.at(kt, 1, i) + Matrex.at(kt, 3, i) <= Matrex.at(it, 1, j)
                   || Matrex.at(it, 1, j) + Matrex.at(it, 3, j) <= Matrex.at(kt, 1, i)
                   || Matrex.at(kt, 2, i) + Matrex.at(kt, 4, i) <= Matrex.at(it, 2, j)
                   || Matrex.at(it, 2, j) + Matrex.at(it, 4, j) <= Matrex.at(kt, 2, i)
              ) do
                vlmt = Matrex.set(vlmt, i, 1, 0)
                {:halt, vlmt}
              else
                {:cont, vlmt}
              end
            end)
          else
            vlmt
          end
        end)
      else
        vlmt
      end
    kt =
      if vlmt do
        Robox.Matrix.get_columns(kt, Robox.Matrix.find(vlmt, 1))
      end
    s =
      if kt do
        Matrex.zeros(kt[:cols], 1)
      end
    xplist = [] #use LIST instead of 3d-matrix, as xplist is a list of matrices
    {xplist, s} =
      if kt do
        Enum.reduce(1..kt[:cols], {xplist, s}, fn i, {xplist, s} ->
          xpWithKt = xp
          xpWithKt = Robox.Gbp2d.Xp.gbp2d_xp_update_xp(bn, it, Robox.Matrix.get_columns(kt, [i]), xpWithKt)
          xplist = List.insert_at(xplist, i, xpWithKt)
          s = Matrex.set(s, i, 1, gbp2d_it_scorer_ktlist(xpWithKt))
          {xplist, s}
        end)
      else
        {xplist, s}
      end
    {kt, xplist, s} = gbp2d_it_purify_ktlist(kt, xplist, s, nlmt)
    %Ktlist2d{
      n: kt_n_cols,
      kt: kt,
      xp: xplist,
      s: s
    }
  end
  
  @spec gbp2d_it_create_ktldht(%Matrex{}) :: struct
  def gbp2d_it_create_ktldht(ktinit) do
    l = Matrex.at(ktinit, 3, 1)
    d = Matrex.at(ktinit, 4, 1)
    vlmt = Matrex.zeros(2, 1)
    ktldht = Matrex.zeros(2, 2)
             |> Matrex.set(1, 1, l)
             |> Matrex.set(2, 1, d)
             |> Matrex.set(1, 2, d)
             |> Matrex.set(2, 2, l)

    vlmt = Matrex.set(vlmt, 1, 1, 1)
    vlmt =
      if l == d do
        vlmt
      else
        Matrex.ones(2, 1)
      end

    Robox.Matrix.get_columns(ktldht, Robox.Matrix.find(vlmt, 1))
  end

  @doc """
  score kt with number of new extreme point and available scales
  """
  @spec gbp2d_it_scorer_ktlist(%Matrex{}) :: number
  def gbp2d_it_scorer_ktlist(xpWithKt) do
    rs = if xpWithKt != nil do
      Enum.reduce(1..xpWithKt[:cols], Matrex.zeros(1, xpWithKt[:cols]), fn i, rs ->
        Matrex.set(rs, 1, i, Matrex.at(xpWithKt, 3, i) * Matrex.at(xpWithKt, 4, i))
      end)
    else
      nil
    end
    rs_sum = if rs != nil do
      Matrex.sum(rs)
    else
      0
    end

    if rs_sum != 0 do
      Enum.reduce(1..xpWithKt[:cols], rs, fn i, rs ->
        Matrex.set(rs, 1, i, Matrex.at(rs, 1, i) / rs_sum)
      end)
      rs_log = Robox.Matrix.log(rs)
      s = Matrex.zeros(1, rs[:cols])
          - Matrex.sum(Enum.reduce(1..rs[:cols], s, fn i, s ->
            rs_elem = Matrex.at(rs, 1, i)
            rs_log_elem = Matrex.at(rs_log, 1, i)
            Matrex.set(s, 1, i, rs_elem * rs_log_elem)
          end))
    else
      0
    end
  end

  @doc """
  purify ktlist with s and nlmt
  """
  @spec gbp2d_it_purify_ktlist(nil, [], nil, number) :: {nil, [], nil}
  def gbp2d_it_purify_ktlist(nil, [], nil, _) do
    {nil, [], nil}
  end

  @spec gbp2d_it_purify_ktlist(%Matrex{}, list, %Matrex{}, number) :: {%Matrex{}, list, %Matrex{}}
  def gbp2d_it_purify_ktlist(kt, xplist, s, nlmt) do
    sidx = Robox.Matrix.sort_index(s) # Enum.with_index in Robox.Matrix.sort_index return indices starts from 0
           |> Matrex.apply(fn i -> i + 1 end) # Matrix indices should start from 1

    if nlmt == 0 || nlmt >= kt[:cols] do
      xplist0 = Enum.map(1..s[:cols], fn i ->
        Enum.at(xplist, trunc(Enum.at(sidx, i - 1) - 1))
      end)

      kt = Robox.Matrix.get_columns(kt, sidx)
      xplist = xplist0
      s = Robox.Matrix.get_rows(s, sidx)

      {kt, xplist, s}
    else
      slmt = Robox.Matrix.linspace(1, nlmt, nlmt)
      xplist0 = Enum.map(1..nlmt, fn i ->
        Enum.at(xplist, trunc(Enum.at(sidx, i - 1) - 1))
      end)
      kt = Robox.Matrix.get_columns(kt, Enum.at(sidx, slmt - 1))
      xplist = xplist0
      s = Robox.Matrix.get_rows(s, Matrex.to_list(Robox.Matrix.get_rows(sidx, slmt)))
      {kt, xplist, s}
    end
  end
end
