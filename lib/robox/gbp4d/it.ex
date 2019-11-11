defmodule Robox.Gbp4d.It do
  alias Robox.{Gbp4d.Ktlist4d, Gbp4d.Xp, Matrix}

  def gbp4d_it_create_ktlist(bn, it, xp, ktinit, nlmt) do
    # init kt ldh
    ktldht = gbp4d_it_create_ktldht(ktinit)

    # init kt
    kt = Matrex.zeros(8, ktldht[:cols] * xp[:cols])

    # ktinit -> kt w.r.t all feasible fit w.r.t orientation and xp position - impose xp residual space as fast filt
    vlmt = Matrex.zeros(ktldht[:cols] * xp[:cols], 1)

    # fast filt w.r.t extreme point residual space
    {_, vlmt, kt} =
      Enum.reduce(1..ktldht[:cols], {1, vlmt, kt}, fn i, {ij, vlmt, kt} ->
        Enum.reduce(1..xp[:cols], {ij, vlmt, kt}, fn j, {ij, vlmt, kt} ->
          if Matrex.at(ktldht, 1, i) <= Matrex.at(xp, 5, j) &&
               Matrex.at(ktldht, 2, i) <= Matrex.at(xp, 6, j) &&
               Matrex.at(ktldht, 3, i) <= Matrex.at(xp, 7, j) &&
               Matrex.at(ktldht, 4, i) <= Matrex.at(xp, 8, j) do
            ij = (i - 1) * xp[:cols] + (j - 1) + 1

            vlmt = Matrex.set(vlmt, ij, 1, 1)

            kt =
              kt
              |> Matrex.set(1, ij, Matrex.at(xp, 1, j))
              |> Matrex.set(2, ij, Matrex.at(xp, 2, j))
              |> Matrex.set(3, ij, Matrex.at(xp, 3, j))
              |> Matrex.set(4, ij, Matrex.at(xp, 4, j))
              |> Matrex.set(5, ij, Matrex.at(ktldht, 1, i))
              |> Matrex.set(6, ij, Matrex.at(ktldht, 2, i))
              |> Matrex.set(7, ij, Matrex.at(ktldht, 3, i))
              |> Matrex.set(8, ij, Matrex.at(ktldht, 4, i))

            {ij, vlmt, kt}
          else
            {ij, vlmt, kt}
          end
        end)
      end)

    kt = Matrix.get_columns(kt, Matrix.find(vlmt, 1))

    # ktinit -> kt w.r.t all feasible fit w.r.t orientation and xp position - impost no conflict with existing it as thorough slow filt
    vlmt =
      if kt do
        Matrex.ones(kt[:cols], 1)
      end

    # slow filt w.r.t no conflict all existing it
    vlmt =
      try do
        if kt == nil do
          throw({:break, vlmt})
        end

        Enum.reduce(1..kt[:cols], vlmt, fn i, vlmt ->
          if it != nil do
            Enum.reduce(1..it[:cols], vlmt, fn j, vlmt ->
              # weight on separate single dimension
              # Matrex.at(kt, 4, i) == sum(Matrex.at(it, 8, j))

              if !(Matrex.at(kt, 1, i) + Matrex.at(kt, 5, i) <= Matrex.at(it, 1, j) ||
                     Matrex.at(it, 1, j) + Matrex.at(it, 5, j) <= Matrex.at(kt, 1, i) ||
                     Matrex.at(kt, 2, i) + Matrex.at(kt, 6, i) <= Matrex.at(it, 2, j) ||
                     Matrex.at(it, 2, j) + Matrex.at(it, 6, j) <= Matrex.at(kt, 2, i) ||
                     Matrex.at(kt, 3, i) + Matrex.at(kt, 7, i) <= Matrex.at(it, 3, j) ||
                     Matrex.at(it, 3, j) + Matrex.at(it, 7, j) <= Matrex.at(kt, 3, i)) do
                throw({:break, Matrex.set(vlmt, i, 1, 0)})
              else
                vlmt
              end
            end)
          else
            vlmt
          end
        end)
      catch
        {:break, vlmt} ->
          vlmt
      end

    kt =
      if vlmt do
        Matrix.get_columns(kt, Matrix.find(vlmt, 1))
      end

    # init ktlist
    s =
      if kt do
        Matrex.zeros(kt[:cols], 1)
      end

    xplist = []

    xp_with_kt = nil

    {_, xplist, s} =
      if kt do
        Enum.reduce(1..kt[:cols], {xp_with_kt, xplist, s}, fn i, {_, xplist, s} ->
          xp_with_kt = Xp.gbp4d_xp_update_xp(bn, it, Matrex.column(kt, i), xp)

          xplist = xplist ++ [xp_with_kt]

          s = Matrex.set(s, i, 1, gbp4d_it_scorer_ktlist(xp_with_kt))

          {xp_with_kt, xplist, s}
        end)
      else
        {xp_with_kt, xplist, s}
      end

    {kt, xplist, s} = gbp4d_it_purify_ktlist(kt, xplist, s, nlmt)

    %Ktlist4d{
      n: if(kt, do: kt[:cols], else: 0),
      kt: kt,
      xp: xplist,
      s: s
    }
  end

  @spec gbp4d_it_purify_ktlist(nil, [], nil, number) :: {nil, [], nil}
  def gbp4d_it_purify_ktlist(nil, [], nil, _) do
    {nil, [], nil}
  end

  @spec gbp4d_it_purify_ktlist(%Matrex{}, list, %Matrex{}, number) :: {%Matrex{}, list, %Matrex{}}
  def gbp4d_it_purify_ktlist(kt, xplist, s, nlmt) do
    sidx = Matrix.sort_index(s)

    if nlmt == 0 || nlmt >= kt[:cols] do
      xplist0 =
        Enum.map(1..s[:rows], fn i ->
          r =
            Matrex.at(sidx, i, 1)
            |> trunc()

          Enum.at(xplist, r)
        end)

      kt = Matrix.get_columns(kt, Matrex.add(sidx, 1))
      xplist = xplist0
      s = Matrix.get_rows(s, Matrex.add(sidx, 1))

      {kt, xplist, s}
    else
      slmt =
        Matrix.linspace(0, nlmt - 1, nlmt)
        |> Matrex.add(1)

      xplist0 =
        Enum.reduce(1..nlmt, [], fn i, xplist0 ->
          xplist0 ++ [Enum.at(xplist, trunc(sidx[i]))]
        end)

      r =
        Matrix.get_rows(sidx, slmt)
        |> Matrex.add(1)

      kt = Matrix.get_columns(kt, r)
      xplist = xplist0

      s = Matrix.get_rows(s, r)

      {kt, xplist, s}
    end
  end

  def gbp4d_it_scorer_ktlist(nil) do
    0.0
  end

  def gbp4d_it_scorer_ktlist(xp_with_kt) do
    # calculcate available space via a single extreme point
    rs =
      Matrex.multiply(Matrex.row(xp_with_kt, 5), Matrex.row(xp_with_kt, 6))
      |> Matrex.multiply(Matrex.row(xp_with_kt, 7))

    # % and / element-wise mutipilication and division
    rs = Matrex.divide(rs, Matrix.sum(rs))

    # score a xpWithKt as the entropy of the available space
    # score s is the smaller the better: prefer the configuration where less and large available dominate
    -Matrix.sum(Matrex.multiply(rs, Matrix.log(rs)))
  end

  def gbp4d_it_create_ktldht(ktinit) do
    # init
    l = ktinit[5]
    d = ktinit[6]
    h = ktinit[7]
    w = ktinit[8]

    # main
    # create all possible none duplicated orientations via mapping:
    # ktinit [l d h w] -> [l d h w; l h d w; d l h w; d h l w; h l d w; h d l w;]
    ktldht =
      Matrex.new([
        [l, l, d, d, h, h],
        [d, h, l, h, l, d],
        [h, d, h, l, d, l],
        [w, w, w, w, w, w]
      ])

    vlmt =
      Matrex.zeros(6, 1)
      |> Matrex.set(1, 1, 1)

    vlmt =
      cond do
        l == d && d == h ->
          vlmt

        l == d ->
          vlmt
          |> Matrex.set(2, 1, 1)
          |> Matrex.set(5, 1, 1)

        d == h ->
          vlmt
          |> Matrex.set(3, 1, 1)
          |> Matrex.set(4, 1, 1)

        h == l ->
          vlmt
          |> Matrex.set(2, 1, 1)
          |> Matrex.set(3, 1, 1)

        true ->
          Matrex.ones(6, 1)
      end

    Matrix.get_columns(ktldht, Matrix.find(vlmt, 1))
  end
end
