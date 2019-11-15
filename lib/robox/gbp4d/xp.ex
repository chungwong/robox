defmodule Robox.Gbp4d.Xp do
  alias Robox.{GbpU, Matrix}

  def gbp4d_xp_create_xp(_bn, _it) do
  end

  @spec gbp4d_xp_update_xp(%Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp4d_xp_update_xp(bn, it, kt, xp) do
    try do
      # init
      xp =
        if it == nil && kt == nil do
          xp =
            Matrex.zeros(8, 1)
            |> Matrex.set(5, 1, bn[1])
            |> Matrex.set(6, 1, bn[2])
            |> Matrex.set(7, 1, bn[3])
            |> Matrex.set(8, 1, bn[4])

          throw({:return, xp})
        else
          xp
        end

      # construct xp input
      # remove extreme points that is taken by or fallen into kt
      # and also update it extreme point residual space w.r.t kt
      xp = gbp4d_xp_update_xp_ikt(it, kt, xp)

      # calculate xp output
      # calculate new extreme points from 6 projections of 3 points
      xp_update = Matrex.fill(8, 6, :nan)

      # maxBound: track 6 projection xp position location
      # from min x-left y-behind z-bottom (w-weight) to max x-right y-front z-top (w-weight)
      # arma::vec maxBound = arma::zeros<arma::vec>(6) - 1; // init maxBound without itBnd save computing cost
      max_bound = Matrex.zeros(6, 1)

      # minBound: track 6 projection xp residual space as
      # from max x-right y-front z-top (w-weight) to min x-left y-behind z-bottom (w-weight)
      min_bound =
        Enum.reduce(1..6, Matrex.zeros(4, 6), fn i, min_bound ->
          min_bound
          # init x-right of all 6 projected extreme point
          |> Matrex.set(1, i, bn[1])
          # init y-front of all 6 projected extreme point
          |> Matrex.set(2, i, bn[2])
          # init z-top   of all 6 projected extreme point
          |> Matrex.set(3, i, bn[3])
          # init w-wlmt  of all 6 projected extreme point
          |> Matrex.set(4, i, bn[4])
        end)

      # calculate xp_update x, y, z, w extreme point position

      {_, xp_update} = gbp4d_xp_update_xp_spg(it, kt, max_bound, xp_update)

      # calculate xpUpdate l, d, h, w residual space along x, y, z, w
      {_, xp_update} = gbp4d_xp_update_rs_spg(it, kt, min_bound, xp_update)

      # prog xp_update remove nan in x, y, z, w, l, d, h, w
      g =
        Enum.reduce(1..7, Matrex.zeros(6, 1), fn i, g ->
          if Matrex.at(xp_update, i, 1) == :nan do
            Matrex.set(g, i, 1, 1)
          else
            g
          end
        end)

      xp_update = Matrix.get_columns(xp_update, Matrix.find(g, 0))

      # join xpUpdate into xp
      # xp is rounded for comparison in unique_cols/1
      xp =
        Matrix.concat(xp, xp_update)
        |> Matrex.apply(fn val ->
          Float.round(val, 2)
        end)
        |> GbpU.unique_cols()
        |> gbp4d_xp_purify_xp()

      # sort xp via non-decreasing order of z, y, x
      ulmt = Matrix.linspace(2, 0, 3)

      GbpU.sort_via_rows(xp, ulmt)
    catch
      {:return, xp} -> xp
    end
  end

  def gbp4d_xp_purify_xp(xp) do
    # remove xp with residual space == 0
    g0 =
      Enum.reduce(1..xp[:cols], Matrex.zeros(xp[:cols], 1), fn i, g0 ->
        if xp[4][i] == 0 || xp[5][i] == 0 || xp[6][i] == 0 || xp[7][i] == 0 do
          Matrex.set(g0, i, 1, 1)
        else
          g0
        end
      end)

    xp = Matrix.get_columns(xp, Matrix.find(g0, 0))

    # remove xp dominated by other xp in list
    # if x, y, z, w < x', y', z', w' and l, d, h, w > l', d', h', w' then
    # (x, y, z, w, l, d, h, w) dominant (x', y', z', w', l', d', h', w')
    # so remove (x', y', z', w', l', d', h', w') from xp list

    if xp do
      g1 =
        Enum.reduce(1..xp[:cols], Matrex.zeros(xp[:cols], 1), fn i, g1 ->
          Enum.reduce(1..xp[:cols], g1, fn j, g1 ->
            # xp[4][i] <= xp[3][j] && // w == w' always true - w and w' are both sum of weight of all it in bn
            # && xp[8][i] >= xp[8][j] // w == w' always true - w and w' are both residual weight available for bn
            if i != j && xp[1][i] <= xp[1][j] && xp[2][i] <= xp[2][j] && xp[3][i] <= xp[3][j] &&
                 xp[5][i] >= xp[5][j] && xp[6][i] >= xp[6][j] && xp[7][i] >= xp[7][j] do
              Matrex.set(g1, j, 1, 1)
            else
              g1
            end
          end)
        end)

      Matrix.get_columns(xp, Matrix.find(g1, 0))
    end
  end

  @spec gbp4d_xp_update_xp_ikt(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp4d_xp_update_xp_ikt(_it, kt, xp) do
    # remove extreme points that is taken by or fallen into kt
    vlmt = Matrex.zeros(xp[:cols], 1)

    vlmt =
      Enum.reduce(1..xp[:cols], vlmt, fn i, vlmt ->
        if kt[1] <= Matrex.at(xp, 1, i) && Matrex.at(xp, 1, i) < kt[1] + kt[5] &&
             kt[2] <= Matrex.at(xp, 2, i) && Matrex.at(xp, 2, i) < kt[2] + kt[6] &&
             kt[3] <= Matrex.at(xp, 3, i) && Matrex.at(xp, 3, i) < kt[3] + kt[7] do
          # kt[4] <= xp[4][i] && xp[4][i] < kt[4] + kt[8] // always true as long as kt(7) > 0
          # kt[4] == xp[4][i] for all i - both weight of all it in bn, kt(7) weight of new kt
          Matrex.set(vlmt, i, 1, 1)
        else
          vlmt
        end
      end)

    xp = Matrix.get_columns(xp, Matrix.find(vlmt, 0))

    if xp != nil do
      # and also update it extreme point residual space w.r.t kt
      Enum.reduce(1..xp[:cols], xp, fn i, xp ->
        xp =
          if Matrex.at(xp, 1, i) <= kt[1] && Matrex.at(xp, 2, i) >= kt[2] &&
               Matrex.at(xp, 2, i) < kt[2] + kt[6] && Matrex.at(xp, 3, i) >= kt[3] &&
               Matrex.at(xp, 3, i) < kt[3] + kt[7] do
            Matrex.set(xp, 5, i, min(Matrex.at(xp, 5, i), kt[1] - Matrex.at(xp, 1, i)))
          else
            xp
          end

        xp =
          if Matrex.at(xp, 2, i) <= kt[2] && Matrex.at(xp, 3, i) >= kt[3] &&
               Matrex.at(xp, 3, i) < kt[3] + kt[7] && Matrex.at(xp, 1, i) >= kt[1] &&
               Matrex.at(xp, 1, i) < kt[1] + kt[5] do
            Matrex.set(xp, 6, i, min(Matrex.at(xp, 6, i), kt[2] - Matrex.at(xp, 2, i)))
          else
            xp
          end

        xp =
          if Matrex.at(xp, 3, i) <= kt[3] && Matrex.at(xp, 1, i) >= kt[1] &&
               Matrex.at(xp, 1, i) < kt[1] + kt[5] && Matrex.at(xp, 2, i) >= kt[2] &&
               Matrex.at(xp, 2, i) < kt[2] + kt[6] do
            Matrex.set(xp, 7, i, min(Matrex.at(xp, 7, i), kt[3] - Matrex.at(xp, 3, i)))
          else
            xp
          end

        xp
        # weight on separate single dimension - weight is holding
        |> Matrex.set(4, i, Matrex.at(xp, 4, i) + kt[8])
        # weight on separate single dimension - weight available
        |> Matrex.set(8, i, Matrex.at(xp, 8, i) - kt[8])
      end)
    else
      xp
    end
  end

  def gbp4d_xp_update_xp_spg(it, kt, max_bound, xp_update) do
    max_bound =
      if it do
        Enum.reduce(1..it[:cols], max_bound, fn i, max_bound ->
          gbp4d_xp_update_maxbnd(Matrex.column(it, i), kt, max_bound)
        end)
      else
        max_bound
      end

    xp_update =
      xp_update
      |> Matrex.set(1, 1, kt[1] + kt[5])
      |> Matrex.set(2, 1, max_bound[1])
      |> Matrex.set(3, 1, kt[3])
      # weight on separate single dimension 
      |> Matrex.set(4, 1, kt[4] + kt[8])
      |> Matrex.set(1, 2, kt[1] + kt[5])
      |> Matrex.set(2, 2, kt[2])
      |> Matrex.set(3, 2, max_bound[2])
      # weight on separate single dimension
      |> Matrex.set(4, 2, kt[4] + kt[8])
      |> Matrex.set(1, 3, kt[1])
      |> Matrex.set(2, 3, kt[2] + kt[6])
      |> Matrex.set(3, 3, max_bound[3])
      # weight on separate single dimension
      |> Matrex.set(4, 3, kt[4] + kt[8])
      |> Matrex.set(1, 4, max_bound[4])
      |> Matrex.set(2, 4, kt[2] + kt[6])
      |> Matrex.set(3, 4, kt[3])
      # weight on separate single dimension
      |> Matrex.set(4, 4, kt[4] + kt[8])
      |> Matrex.set(1, 5, max_bound[5])
      |> Matrex.set(2, 5, kt[2])
      |> Matrex.set(3, 5, kt[3] + kt[7])
      # weight on separate single dimension
      |> Matrex.set(4, 5, kt[4] + kt[8])
      |> Matrex.set(1, 6, kt[1])
      |> Matrex.set(2, 6, max_bound[6])
      |> Matrex.set(3, 6, kt[3] + kt[7])
      # weight on separate single dimension
      |> Matrex.set(4, 6, kt[4] + kt[8])

    {max_bound, xp_update}
  end

  def gbp4d_xp_update_rs_spg(it, kt, min_bound, xp_update) do
    min_bound =
      if it do
        Enum.reduce(1..it[:cols], min_bound, fn i, min_bound ->
          gbp4d_xp_update_minbnd(Matrex.column(it, i), kt, min_bound, xp_update)
        end)
      else
        min_bound
      end

    xp_update =
      Enum.reduce(1..6, xp_update, fn i, xp_update ->
        xp_update
        |> Matrex.set(5, i, min_bound[1][i] - xp_update[1][i])
        |> Matrex.set(6, i, min_bound[2][i] - xp_update[2][i])
        |> Matrex.set(7, i, min_bound[3][i] - xp_update[3][i])
        # weight on separate single dimension
        |> Matrex.set(8, i, min_bound[4][i] - xp_update[4][i])
      end)

    {min_bound, xp_update}
  end

  def gbp4d_xp_update_minbnd(it, _kt, min_bound, xp_update) do
    # projecting kt xp -> it on reverse direction w.r.t gbp3d_xp_it_pjt_kt.png
    # arma::uvec ik = gbp4d_xp_it_qjt_kt(it, kt);

    # construct virtual kt as xpUpdate with l = 0, d = 0, h = 0, w = 0
    # since residual space of extreme point is related to the x, y, z, w of exteme point itself
    akt = Matrex.zeros(8, 1)
    aik = Matrex.zeros(6, 1)

    {_, _, min_bound} =
      Enum.reduce(1..6, {akt, aik, min_bound}, fn i, {akt, _aik, min_bound} ->
        # init
        # init residual space without creating another itBnd and save computation cost (skip)
        # xpUpdate(4, i) = minBound(0, i) - xpUpdate(0, i);
        # xpUpdate(5, i) = minBound(1, i) - xpUpdate(1, i);
        # xpUpdate(6, i) = minBound(2, i) - xpUpdate(2, i);
        # xpUpdate(7, i) = minBound(3, i) - xpUpdate(3, i);

        # create a virtual kt as a single point with l = 0, d = 0, h = 0, w = 0
        akt =
          akt
          |> Matrex.set(1, 1, xp_update[1][i])
          |> Matrex.set(2, 1, xp_update[2][i])
          |> Matrex.set(3, 1, xp_update[3][i])
          |> Matrex.set(4, 1, xp_update[4][i])
          |> Matrex.set(5, 1, 0.00)
          |> Matrex.set(6, 1, 0.00)
          |> Matrex.set(7, 1, 0.00)
          |> Matrex.set(8, 1, 0.00)

        # projecting kt xp -> it on reverse direction w.r.t gbp3d_xp_it_pjt_kt.png
        aik = gbp4d_xp_it_qjt_kt(it, akt)

        # since akt(4) = 0.00; akt(5) = 0.00; akt(6) = 0.00; =>
        # aik(0) == aik(5); aik(1) == aik(2); aik(3) == aik(4);

        # it block on the way from extreme point to x-right
        min_bound =
          if aik[4] > 0 && aik[5] > 0 do
            Matrex.set(min_bound, 1, i, min(it[1], min_bound[1][i]))
          else
            min_bound
          end

        min_bound =
          if aik[6] > 0 && aik[1] > 0 do
            Matrex.set(min_bound, 2, i, min(it[2], min_bound[2][i]))
          else
            min_bound
          end

        min_bound =
          if aik[2] > 0 && aik[3] > 0 do
            Matrex.set(min_bound, 3, i, min(it[3], min_bound[3][i]))
          else
            min_bound
          end

        {akt, aik, min_bound}
      end)

    min_bound
  end

  def gbp4d_xp_update_maxbnd(it, kt, max_bound) do
    # projecting kt xp -> it along with direction w.r.t gbp4d_xp_it_pjt_kt.png
    ik = gbp4d_xp_it_pjt_kt(it, kt)

    # direction x-Y: kt-x-corner-move project-toward->Y
    max_bound =
      if ik[1] > 0 && it[2] + it[6] > max_bound[1] do
        Matrex.set(max_bound, 1, 1, it[2] + it[6])
      else
        max_bound
      end

    # direction x-Z: kt-x-corner-move project-toward->Z
    max_bound =
      if ik[2] > 0 && it[3] + it[7] > max_bound[2] do
        Matrex.set(max_bound, 2, 1, it[3] + it[7])
      else
        max_bound
      end

    # direction y-Z: kt-y-corner-move project-toward->Z
    max_bound =
      if ik[3] > 0 && it[3] + it[7] > max_bound[3] do
        Matrex.set(max_bound, 3, 1, it[3] + it[7])
      else
        max_bound
      end

    # direction y-X: kt-y-corner-move project-toward-X
    max_bound =
      if ik[4] > 0 && it[1] + it[5] > max_bound[4] do
        Matrex.set(max_bound, 4, 1, it[1] + it[5])
      else
        max_bound
      end

    # direction z-X: kt-z-corner-move project-toward->X
    max_bound =
      if ik[5] > 0 && it[1] + it[5] > max_bound[5] do
        Matrex.set(max_bound, 5, 1, it[1] + it[5])
      else
        max_bound
      end

    # direction z-Y: kt-z-corner-move project-toward->Y
    if ik[6] > 0 && it[2] + it[6] > max_bound[6] do
      Matrex.set(max_bound, 6, 1, it[2] + it[6])
    else
      max_bound
    end
  end

  def gbp4d_xp_it_qjt_kt(it, kt) do
    ik = []

    # direction x-Y
    ik =
      ik ++
        [
          kt[2] + kt[6] <= it[2] && it[1] <= kt[1] + kt[5] && kt[1] + kt[5] < it[1] + it[5] &&
            it[3] <= kt[3] && kt[3] < it[3] + it[7]
        ]

    # direction x-Z
    ik =
      ik ++
        [
          kt[3] + kt[7] <= it[3] && it[1] <= kt[1] + kt[5] && kt[1] + kt[5] < it[1] + it[5] &&
            it[2] <= kt[2] && kt[2] < it[2] + it[6]
        ]

    # direction y-Z
    ik =
      ik ++
        [
          kt[3] + kt[7] <= it[3] && it[2] <= kt[2] + kt[6] && kt[2] + kt[6] < it[2] + it[6] &&
            it[1] <= kt[1] && kt[1] < it[1] + it[5]
        ]

    # direction y-X
    ik =
      ik ++
        [
          kt[1] + kt[5] <= it[1] && it[2] <= kt[2] + kt[6] && kt[2] + kt[6] < it[2] + it[6] &&
            it[3] <= kt[3] && kt[3] < it[3] + it[7]
        ]

    # direction z-X
    ik =
      ik ++
        [
          kt[1] + kt[5] <= it[1] && it[3] <= kt[3] + kt[7] && kt[3] + kt[7] < it[3] + it[7] &&
            it[2] <= kt[2] && kt[2] < it[2] + it[6]
        ]

    # direction z-Y
    ik =
      ik ++
        [
          kt[2] + kt[6] <= it[2] && it[3] <= kt[3] + kt[7] && kt[3] + kt[7] < it[3] + it[7] &&
            it[1] <= kt[1] && kt[1] < it[1] + it[5]
        ]

    ik =
      Enum.map(ik, fn v ->
        [Matrix.boolean_to_integer(v)]
      end)

    Matrex.new(ik)
  end

  def gbp4d_xp_it_pjt_kt(it, kt) do
    ik = []

    # direction x-Y
    ik =
      ik ++
        [
          it[2] + it[6] <= kt[2] && it[1] <= kt[1] + kt[5] && kt[1] + kt[5] < it[1] + it[5] &&
            it[3] <= kt[3] && kt[3] < it[3] + it[7]
        ]

    # direction x-Z
    ik =
      ik ++
        [
          it[3] + it[7] <= kt[3] && it[1] <= kt[1] + kt[5] && kt[1] + kt[5] < it[1] + it[5] &&
            it[2] <= kt[2] && kt[2] < it[2] + it[6]
        ]

    # direction y-Z
    ik =
      ik ++
        [
          it[3] + it[7] <= kt[3] && it[2] <= kt[2] + kt[6] && kt[2] + kt[6] < it[2] + it[6] &&
            it[1] <= kt[1] && kt[1] < it[1] + it[5]
        ]

    # direction y-X
    ik =
      ik ++
        [
          it[1] + it[5] <= kt[1] && it[2] <= kt[2] + kt[6] && kt[2] + kt[6] < it[2] + it[6] &&
            it[3] <= kt[3] && kt[3] < it[3] + it[7]
        ]

    # direction z-X
    ik =
      ik ++
        [
          it[1] + it[5] <= kt[1] && it[3] <= kt[3] + kt[7] && kt[3] + kt[7] < it[3] + it[7] &&
            it[2] <= kt[2] && kt[2] < it[2] + it[6]
        ]

    # direction z-Y
    ik =
      ik ++
        [
          it[2] + it[6] <= kt[2] && it[3] <= kt[3] + kt[7] && kt[3] + kt[7] < it[3] + it[7] &&
            it[1] <= kt[1] && kt[1] < it[1] + it[5]
        ]

    ik =
      Enum.map(ik, fn v ->
        [Matrix.boolean_to_integer(v)]
      end)

    Matrex.new(ik)
  end
end
