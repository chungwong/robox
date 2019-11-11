defmodule Robox.Gbp2d.Xp do
  @moduledoc """
  """

  @doc """
  """
  @spec gbp2d_xp_update_xp(%Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_update_xp(bn, it, kt, xp) do
    if it == nil && kt == nil do
      xp = Matrex.zeros(4, 1)
      xp = Matrex.set(xp, 3, 1, Matrex.at(bn, 1, 1))
      Matrex.set(xp, 4, 1, Matrex.at(bn, 2, 1))
    else
      xp = gbp2d_xp_update_xp_ikt(kt, xp);
      xpUpdate = Matrex.fill(4, 2, :nan)
      maxBound = Matrex.zeros(2, 1)
      minBound = Enum.reduce(1..2, Matrex.zeros(2, 2), fn i, minBound ->
        minBound = Matrex.set(minBound, 1, i, Matrex.at(bn, 1, 1))
        Matrex.set(minBound, 2, i, Matrex.at(bn, 2, 1))
      end)
      # calculate xpUpdate x, y extreme point position
      {xpUpdate, _} = gbp2d_xp_update_xp_spg(it, kt, maxBound, xpUpdate)
      # calculate xpUpdate l, d residual space along x, y
      {xpUpdate, _} = gbp2d_xp_update_rs_spg(it, minBound, xpUpdate)
      g = Matrex.zeros(2, 1)
      g = Enum.reduce(1..2, g, fn i, g ->
            if Matrex.contains?(Matrex.row(g, i), :nan) do
              Matrex.set(g, i, 1, 1)
            else
              g
            end
          end)
      xpUpdate = Robox.Matrix.get_columns(xpUpdate, Robox.Matrix.find(g, 0))
      xp =
        if xp == nil do # Matrex.concat cannot handle nil
          Robox.GbpU.unique_cols(xpUpdate)
        else
          Robox.GbpU.unique_cols(Matrex.concat(xp, xpUpdate))
        end
      xp = gbp2d_xp_purify_xp(xp)
      # sort xp via non-decreasing order of y, x
      ulmt = Matrex.new("1;0")
      xp = Robox.GbpU.sort_via_rows(xp, ulmt)
      xp
    end
  end

  @doc """
  update current extreme point xp list based on it w.r.t fit kt into bn
  """
  @spec gbp2d_xp_update_xp_ikt(%Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_update_xp_ikt(kt, xp) do
    try do
      if xp == nil do
        throw(:return)
      end

      vlmt = Matrex.zeros(xp[:columns], 1)
      vlmt = Enum.reduce(1..xp[:columns], vlmt, fn i, vlmt ->
        if(Matrex.at(kt, 1, 1) <= Matrex.at(xp, 1, i)
          && Matrex.at(xp, 1, i) < Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1)
          && Matrex.at(kt, 1, 1) <= Matrex.at(xp, 1, i)
          && Matrex.at(xp, 1, i) < Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1)
        ) do
          Matrex.set(vlmt, i, 1, 1)
        else
          vlmt
        end
      end)

      xp = Robox.Matrix.get_columns(xp, Robox.Matrix.find(vlmt, 0))

      if xp == nil do
        throw(:return)
      end
      Enum.reduce(1..xp[:columns], xp, fn i, xp ->
        xp = if(Matrex.at(xp, 1, i) <= Matrex.at(kt, 1, 1) &&
                Matrex.at(xp, 2, i) <= Matrex.at(kt, 2, 1) &&
                Matrex.at(xp, 2, i) <= Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1)
        ) do
          Matrex.set(xp, 3, i, Enum.min([Matrex.at(xp, 4, i), Matrex.at(kt, 2, 1) - Matrex.at(xp, 2, i)]))
        else
          xp
        end

        if(Matrex.at(xp, 2, i) <= Matrex.at(kt, 2, 1) &&
                Matrex.at(xp, 1, i) <= Matrex.at(kt, 1, 1) &&
                Matrex.at(xp, 1, i) <= Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1)
        ) do
          Matrex.set(xp, 4, i, min(Matrex.at(xp, 4, i), Matrex.at(kt, 2, 1) - Matrex.at(xp, 2, i)))
        else
          xp
        end
      end)
    catch
      :return -> nil
    end
  end

  @doc """
  calculate xpUpdate x, y extreme point position
  """
  @spec gbp2d_xp_update_xp_spg(%Matrex{}, %Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_update_xp_spg(it, kt, maxBound, xpUpdate) do
    maxBound =
      if it != nil do
        Enum.reduce(1..it[:columns], maxBound, fn i, maxBound ->
          gbp2d_xp_update_maxbnd(Robox.Matrix.get_columns(it, [i]), kt, maxBound)
        end)
      else
        maxBound
      end
    xpUpdate =
      Matrex.set(xpUpdate, 1, 1, Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1))
      |> Matrex.set(2, 1, Matrex.at(maxBound, 1, 1))
      |> Matrex.set(1, 2, Matrex.at(maxBound, 2, 1))
      |> Matrex.set(2, 2, Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1))
    {xpUpdate, maxBound}
  end

  @doc """
  calculate residual space of projected kt xp over each single it in bin
  """
  @spec gbp2d_xp_update_rs_spg(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_update_rs_spg(it, minBound, xpUpdate) do
    minBound = if it != nil do
      Enum.reduce(1..it[:columns], minBound, fn i, minBound ->
        gbp2d_xp_update_minbnd(Robox.Matrix.get_columns(it, [i]), minBound, xpUpdate)
      end)
    else
      minBound
    end
    xpUpdate = Enum.reduce(1..2, xpUpdate, fn i, xpUpdate ->
      xpUpdate = Matrex.set(xpUpdate, 3, i, Matrex.at(minBound, 1, i) - Matrex.at(xpUpdate, 1, i))
      Matrex.set(xpUpdate, 4, i, Matrex.at(minBound, 2, i) - Matrex.at(xpUpdate, 2, i))
    end)
    {xpUpdate, minBound}
  end

  @doc """
  calculate residual space of projected kt xp over each single it in bin
  """
  @spec gbp2d_xp_update_minbnd(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_update_minbnd(it, minBound, xpUpdate) do
    akt = Matrex.zeros(4, 1)

    Enum.reduce(1..2, minBound, fn i, minBound ->
      akt = Matrex.set(akt, 1, 1, Matrex.at(xpUpdate, 1, i))
      akt = Matrex.set(akt, 2, 1, Matrex.at(xpUpdate, 2, i))
      akt = Matrex.set(akt, 3, 1, 0)
      akt = Matrex.set(akt, 4, 1, 0)
      aik = gbp2d_xp_it_qjt_kt(it, akt)
      minBound =
        if Matrex.at(aik, 2, 1) > 0 do
          Matrex.set(minBound, 1, i, Enum.min([Matrex.at(it, 1, 1), Matrex.at(minBound, 1, i)]))
        else
          minBound
        end

        if Matrex.at(aik, 1, 1) > 0 do
          Matrex.set(minBound, 2, i, Enum.min([Matrex.at(it, 2, 1), Matrex.at(minBound, 2, i)]))
        else
          minBound
        end
    end)
  end

  @doc """
  """
  @spec gbp2d_xp_it_qjt_kt(%Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_it_qjt_kt(it, kt) do
    ik = Matrex.zeros(2, 1)

    ik = Matrex.set(ik, 1, 1, boolean_to_number(
         Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1) <= Matrex.at(it, 2, 1)
      && Matrex.at(it, 1, 1)                       <= Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1)
      && Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1) <= Matrex.at(it, 1, 1) + Matrex.at(it, 3, 1)
    ))

    ik = Matrex.set(ik, 2, 1, boolean_to_number(
         Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1) <= Matrex.at(it, 1, 1)
      && Matrex.at(it, 2, 1)                       <= Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1)
      && Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1) <= Matrex.at(it, 2, 1) + Matrex.at(it, 4, 1)
    ))
    ik
  end

  @doc """
  # pure xp list via remove xp with residual space == 0
  # and via remove xp dominated by other xp in the list
  """
  @spec gbp2d_xp_purify_xp(%Matrex{}) :: %Matrex{}
  def gbp2d_xp_purify_xp(xp) do
    xp =
    try do
      if xp == nil do
        throw(:return)
      else
        g0 = Matrex.zeros(xp[:cols], 1)
        g0 = Enum.reduce(1..xp[:cols], g0, fn i, g0 ->
          if Matrex.at(xp, 3, i) == 0 || Matrex.at(xp, 4, i) == 0 do
            Matrex.set(g0, i, 1, 1)
          else
            g0
          end
        end)
        xp = Robox.Matrix.get_columns(xp, Robox.Matrix.find(g0, 0))
        if xp != nil do
          g1 = Matrex.zeros(xp[:cols], 1)
          Enum.reduce(1..xp[:cols], g1, fn i, g1 ->
            Enum.reduce(1..xp[:cols], g1, fn j, g1 ->
              if i != j &&
                Matrex.at(xp, 1, i) <= Matrex.at(xp, 1, j) &&
                Matrex.at(xp, 2, i) <= Matrex.at(xp, 2, j) &&
                Matrex.at(xp, 3, i) <= Matrex.at(xp, 3, j) &&
                Matrex.at(xp, 4, i) <= Matrex.at(xp, 4, j) do
                  Matrex.set(g1, j, 1, 1)
              else
                g1
              end
            end)
          end)
          Robox.Matrix.get_columns(xp, Robox.Matrix.find(g1, 0))
        end
      end
    catch
      :return -> nil
    end
    xp
  end

  @doc """
  """
  @spec gbp2d_xp_update_maxbnd(%Matrex{}, %Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_update_maxbnd(it, kt, maxBound) do
    ik = gbp2d_xp_it_pjt_kt(it, kt)
    maxBound = if Matrex.at(ik, 1, 1) > 0 && Matrex.at(it, 2, 1) + Matrex.at(it, 4, 1) > Matrex.at(maxBound, 1, 1) do
      Matrex.set(maxBound, 1, 1, Matrex.at(it, 2, 1) + Matrex.at(it, 4, 1))
    else
      maxBound
    end
    maxBound = if Matrex.at(ik, 2, 1) && Matrex.at(it, 1, 1) + Matrex.at(it, 3, 1) > Matrex.at(maxBound, 2, 1) do
      Matrex.set(maxBound, 2, 1, Matrex.at(it, 1, 1) + Matrex.at(it, 3, 1))
    else
      maxBound
    end
    maxBound
  end

  @doc """
  """
  @spec boolean_to_number(boolean) :: number
  def boolean_to_number(bool) do
    if(bool == true) do
      1
    else
      0
    end
  end

  @doc """
  can item it take projection of new extreme point xp created by next item kt on the one direction w.r.t gbp2d_xp_it_pjt_kt.png
  """
  @spec gbp2d_xp_it_pjt_kt(%Matrex{}, %Matrex{}) :: %Matrex{}
  def gbp2d_xp_it_pjt_kt(it, kt) do
    ik = Matrex.zeros(2,1)

    ik = Matrex.set(ik, 1, 1, boolean_to_number(
         Matrex.at(it, 2, 1) + Matrex.at(it, 4, 1) <= Matrex.at(kt, 2, 1)
      && Matrex.at(it, 1, 1)                       <= Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1)
      && Matrex.at(kt, 1, 1) + Matrex.at(kt, 3, 1) <= Matrex.at(it, 1, 1) + Matrex.at(it, 3, 1))
    )
    
    Matrex.set(ik, 2, 1, boolean_to_number(
         Matrex.at(it, 1, 1) + Matrex.at(it, 3, 1) <= Matrex.at(kt, 1, 1)
      && Matrex.at(it, 2, 1)                       <= Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1)
      && Matrex.at(kt, 2, 1) + Matrex.at(kt, 4, 1) <= Matrex.at(it, 2, 1) + Matrex.at(it, 4, 1))
    )
  end
end
